function [Theta_pdx,libtree] = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,m_x,m_t,sub_inds,dim,scales,xs,customf,customconv)

    dx = mean(cellfun(@(x)mean(diff(x)),xs(1:end-1)));
    dt = mean(diff(xs{end}));
    Ns = size(U_obs{1});
    if isempty(scales)
        scales = ones(n+dim,1);
    end

    Cfs_ffts = cell(dim,1);
    [mm,nn] = size(Cfs_x);
    for k=1:dim-1
        Cfs_ffts{k} = [zeros(mm,Ns(k)-nn) (m_x*dx*scales(n+k)).^(-(0:mm-1)').*Cfs_x];
        Cfs_ffts{k} = fft(Cfs_ffts{k},[],2);
    end
    [mm,nn] = size(Cfs_t);
    Cfs_ffts{dim} = [zeros(mm,Ns(dim)-nn) (m_t*dt*scales(n+dim)).^(-(0:mm-1)').*Cfs_t];
    Cfs_ffts{dim} = fft(Cfs_ffts{dim},[],2);

    libtree = liblisttree(lib_list(1:end-length(customf)-length(customconv),:),n);

    Theta_pdx = zeros(prod(cellfun(@(x)length(x),sub_inds)),size(lib_list,1));

    for i=1:length(libtree)
        fcn = ones(size(U_obs{1}));
        tags_root = zeros(1,n);
        indout = libtree{i};
        while ~isempty(indout)
            tags = lib_list(indout(1),1:n);
            sametags = indout(ismember(lib_list(indout,1:n),tags,'rows'));
            for k=1:n
                if isreal(tags(k))
                    fcn = fcn.*(U_obs{k} / scales(k)).^(tags(k)-tags_root(k));
                else
                    if imag(tags(k))<0
                        fcn = sin(abs(imag(tags(k)))*U_obs{k});
                    else
                        fcn = cos(abs(imag(tags(k)))*U_obs{k});
                    end
                end
            end
            for ind = sametags
                test_conv_cell = {};
                for k=1:dim
                    test_conv_cell{k} = Cfs_ffts{k}(lib_list(ind,n+k)+1,:);
                end
                
                fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,2);
                Theta_pdx(:,ind) = fcn_conv(:);
            end
            indout = indout(~ismember(indout,sametags));
            tags_root = tags;
        end
    end
    
% ind =1;
% while ind<=size(lib_list,1)-length(customf)-length(customconv)
%     tags = lib_list(ind,1:n);
%     fcn = ones(size(U_obs{1}));
%     for k=1:n
%         if isreal(tags(k))
%             fcn = fcn.*((U_obs{k} / scales(k)).^tags(k));
%         else
%             if imag(tags(k))<0
%                 fcn = sin(abs(imag(tags(k)))*U_obs{k});
%             else
%                 fcn = cos(abs(imag(tags(k)))*U_obs{k});
%             end
%         end
%     end
%     while all(lib_list(ind,1:n) == tags)
%         test_conv_cell = {};
%         for k=1:dim
%             test_conv_cell{k} = Cfs_ffts{k}(lib_list(ind,n+k)+1,:);
%         end
%         fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,2);
%         Theta_pdx(:,ind) = fcn_conv(:);
%         ind = ind+1;
%         if ind > size(lib_list,1)
%             break
%         end
%     end
% end

    for i=1:length(xs)
        L = length(xs{i});
        empty_inds = repmat({1},1,length(xs));
        empty_inds{i} = L;
        xs{i} = reshape(xs{i},empty_inds{:});
    end

    Theta_customf = zeros(size(Theta_pdx,1),length(customf));
    for ind =1:length(customf)    
        u_tags = customf{ind}{1};
        x_tags = customf{ind}{2};
        psi_tags = customf{ind}{3};

        fcn = (U_obs{1} / scales(1)).^u_tags(1);
        for k=2:n
            fcn = fcn.*((U_obs{k} / scales(k)).^u_tags(k));
        end
        for k=1:dim
            fcn = fcn.*((xs{k} / scales(n+k)).^x_tags(k));
        end

        test_conv_cell = {};
        for k=1:dim
            test_conv_cell{k} = Cfs_ffts{k}(psi_tags(k)+1,:);
        end
        fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,2);
        Theta_customf(:,ind) = fcn_conv(:);
    end
    Theta_pdx(:,1+size(lib_list,1)-length(customf)-length(customconv):size(lib_list,1)-length(customconv)) = Theta_customf;

    Theta_conv = zeros(size(Theta_pdx,1),length(customconv));
    for ind=1:length(customconv)
        u_tags = customconv{ind}{1};
        K = customconv{ind}{2};
        psi_tags = customconv{ind}{3};
        v_tags = customconv{ind}{4};
        svdtol = customconv{ind}{5};
        fcn = (U_obs{1} / scales(1)).^u_tags(1);
        for k=2:n
            if u_tags(k)~=0
                fcn = fcn.*((U_obs{k} / scales(k)).^u_tags(k));
            end
        end
        fcn = convNDfftsketch(fcn,xs,K,svdtol);
        for k=1:n
            if v_tags(k)~=0
                fcn = fcn.*((U_obs{k} / scales(k)).^v_tags(k));
            end
        end
        test_conv_cell = {};
        for k=1:dim
            test_conv_cell{k} = Cfs_ffts{k}(psi_tags(k)+1,:);
        end
        fcn = convNDfft(fcn,test_conv_cell,sub_inds,2);
        Theta_conv(:,ind) = fcn(:);
    end
    Theta_pdx(:,1+size(lib_list,1)-length(customconv):end) = Theta_conv;
    % if and(dim==2,length(customconv)>0)
    %     m = [1 1 zeros(1,length(customconv)-2)];
    %     for j=1:length(customconv)/2-1;
    %         m = [m;circshift(m(end,:),2)];
    %     end
    %     m = m';
    %     Theta_pdx = [Theta_pdx(:,1:end-length(customconv)) Theta_pdx(:,end-length(customconv)+1:end)*m];
    % end

end
