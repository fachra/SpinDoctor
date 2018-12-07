function [TOUT,YOUT,MT,difftime,elapsed_time] ...
    = solve_magnetization(gdir,qvalues,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol,...
    Ncmpt,Nboundary,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts,...
    Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,kappa_vec,IC_cmpts,ic_llimit,ic_ulimit)

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG 
global BDELTA SDELTA SEQ OGSEPER 
  
yes = 1;  no = 0;
ODEsolve_atol = atol;
ODEsolve_rtol = rtol;




UG = gdir';
UG = UG/norm(UG);

if(Ncmpt ==1 | abs(max(kappa_vec)) <= 1e-16)
    DO_COUPLING = no;
else
    DO_COUPLING = yes;
end


nexperi = length(sdeltavec);

disp(['Simulating ',num2str(nexperi),' experiments']);

nb = size(qvalues,2);

elapsed_time=zeros(nb, nexperi);

if (1 == 1)
    for icmpt = 1:Ncmpt
        %disp(['Working on cmpt ', num2str(icmpt)]);
        coordinates = Pts_cmpt_reorder{icmpt};
        elements = Ele_cmpt_reorder{icmpt};
        facets = [];
        GX = -sqrt(-1)*UG*coordinates;
        FEM_MAT{icmpt}.Q = sparse(length(coordinates),length(coordinates));
        for iboundary = 1:Nboundary
            neumann = Fac_boundary_reorder{icmpt}{iboundary}';
            neumann_nodes=unique(neumann);
            coeffs_flux_matrix=zeros(max(neumann_nodes),1);
            coeffs_flux_matrix(neumann_nodes)=kappa_vec(iboundary);
            if ~isempty(neumann)
                FEM_MAT{icmpt}.Q = FEM_MAT{icmpt}.Q + flux_matrixP1_3D(neumann,coordinates',coeffs_flux_matrix);
            end;
        end;
        [FEM_MAT{icmpt}.K,volumes]=stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
        FEM_MAT{icmpt}.M=mass_matrixP1_3D(elements',volumes);
        FEM_MAT{icmpt}.A=mass_matrixP1_3D(elements',volumes,GX');
    end
else   
    for icmpt = 1:Ncmpt
        figure;
        model_orig{icmpt} = createpde();
        geometryFromMesh(model_orig{icmpt},Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
        
        %pdeplot3D(model_orig{icmpt}.Mesh,'FaceAlpha',0.1);
        %generateMesh(model_orig{icmpt},'GeometricOrder','linear');
        specifyCoefficients(model_orig{icmpt},'m',0,'d',1,'c',DIFF_cmpts(icmpt),'a',@BTfreqterm_notime_nob,'f',0);
        nface = model_orig{icmpt}.Geometry.NumFaces; %number of face
        %     normal dot (c grad u) + q*u = g
        for iface = 1:nface
            applyBoundaryCondition(model_orig{icmpt},'neumann','face',iface,'g',0,'q',kappa_vec(1),'Vectorized','on'); % 3-D geometry
        end
        FEM_MAT{icmpt} = assembleFEMatrices(model_orig{icmpt});
    end
  
end
for icmpt=1:Ncmpt
    % calculate the IC
    IC_Pts{icmpt} = IC_cmpts(icmpt)*ones(size(Pts_cmpt_reorder{icmpt},2),1);
    ii = find(Pts_cmpt_reorder{icmpt}(1,:)<ic_llimit(1,1) | Pts_cmpt_reorder{icmpt}(1,:)>ic_ulimit(1,1) | ...
        Pts_cmpt_reorder{icmpt}(2,:)<ic_llimit(2,1) | Pts_cmpt_reorder{icmpt}(2,:)>ic_ulimit(2,1) | ...
        Pts_cmpt_reorder{icmpt}(3,:)<ic_llimit(3,1) | Pts_cmpt_reorder{icmpt}(3,:)>ic_ulimit(3,1));
    IC_Pts{icmpt}(ii,1) = 0;
end

if (DO_COUPLING == yes)
    tic
    [FEMcouple_MAT,FEMcouple_ind0,FEMcouple_indf] ...
        = generate_FEM_coupling(FEM_MAT,Ncmpt,Nboundary,...
        Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder);
    %stop
    toc
    % IC for coupled case
    IC_couple=zeros(size(FEMcouple_MAT.M,1),1);
    for icmpt=1:Ncmpt
        IC_couple(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),1) = IC_Pts{icmpt};
    end
    
else
    FEMcouple_MAT = [];
    FEMcouple_ind0 = [];
    FEMcouple_indf = [];    
end

%% solve ODE
for iexperi = 1:nexperi
    
    SDELTA = sdeltavec(iexperi);
    BDELTA = bdeltavec(iexperi);
    TE = SDELTA+BDELTA;
    SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
    omega = 2*pi*npervec(iexperi)/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE
    
    disp(['EXPERIMENT ',num2str(iexperi)]);
    disp([SDELTA,BDELTA,npervec(iexperi)]);
    
    TLIST = [0,SDELTA+BDELTA];
    for ib = 1:nb
        b_start_time = clock;
        % global variable setting QVAL for ODE time stepping
        QVAL = qvalues(iexperi,ib);        
        disp(['Q-VALUE:',mynum2str(QVAL)]);
        
        difftime(iexperi) = seqdifftime;
        
        %% Solving for case of coupling between compartments.
        
        if (DO_COUPLING == yes)
            
            FEM_M = FEMcouple_MAT.M;
            FEM_K = FEMcouple_MAT.K;
            FEM_A = FEMcouple_MAT.A; %*seqprofile(state.time)*QVAL;
            FEM_Q = FEMcouple_MAT.Q;
            
            FEM_G = sparse(zeros(size(FEM_M,1),1));
            
            options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                'Jacobian',@odejac_bt_includeb);            
            disp('***Coupled: start ode solve ode23t'); tic
            sol= ode23t(@odefun_bt_includeb,TLIST,IC_couple,options);
            disp('***Coupled: end ode solve ode23t'); toc
            for icmpt = 1:Ncmpt
                YOUT{iexperi}{ib}{icmpt} = sol.y(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),:);
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
            end
        else
            %% Solving for case of no coupling between compartments.            
            for icmpt = 1:Ncmpt
                FEM_M = FEM_MAT{icmpt}.M;
                FEM_K = FEM_MAT{icmpt}.K;
                FEM_A = FEM_MAT{icmpt}.A;
                FEM_Q = FEM_MAT{icmpt}.Q;
                FEM_G = sparse(zeros(size(FEM_M,1),1));
                options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                    'Jacobian',@odejac_bt_includeb);
                disp('***Uncoupled: start ode solver ode23t'); tic
                
                ICC = IC_Pts{icmpt};
                if (max(abs(ICC))<=1e-16)
                    sol.y = zeros(size(ICC,1),size(TLIST,2));
                    sol.x = TLIST;
                else
                    sol = ode23t(@odefun_bt_includeb,TLIST,ICC,options);
                end
                disp('***Uncoupled: end ode solver ode23t'); toc
                YOUT{iexperi}{ib}{icmpt} = sol.y;
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
            end            
        end
        elapsed_time(ib, iexperi)=etime(clock, b_start_time);
    end
end





