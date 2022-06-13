function true_nz_weights = tags2axi(true_nz_weight_tags,lib_list)
   if ~isempty(true_nz_weight_tags)
       if all(cellfun(@(x) ~isempty(x), true_nz_weight_tags))
          neq = length(true_nz_weight_tags);
          m = size(lib_list,1);
          true_nz_weights = zeros(m,neq);
          for k=1:neq
              [~,loc] = ismember(true_nz_weight_tags{k}(:,1:end-1),lib_list,'rows');
              for j=find(loc)
                  true_nz_weights(loc(j),k) = true_nz_weight_tags{k}(j,end);
              end
              if any(loc==0)
                  disp('true terms missing from lib')
              end
          end
       else
          true_nz_weights=[];
       end
   else
       true_nz_weights=[];
   end
end