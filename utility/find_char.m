function ind=find_char(wholeset,subset,UniformOutput,cmp,verbose)
%%
%     Author: Ying Wang, Min Li
%     Create Time: May 22, 2021
%     Copyright(c): 2020-2022 Ying Wang, yingwangrigel@gmail.com, Min Li, minli.231314@gmail.com
%     Joint China-Cuba LAB, UESTC

if nargin<3 || isempty(UniformOutput)
    UniformOutput=true;
end

if ischar(subset)
    subset={subset};
end


% if nargin<4
%     ignorecase=false;
% end
% if ignorecase
%     cmp=@strcmpi;
% else
%     cmp=@strcmp;
% end
if nargin<4 || isempty(cmp)
    cmp=@strcmp;
end
if nargin<5 || isempty(verbose)
    verbose=false;
end


if UniformOutput
    try
        %         index=cellfun(@(x) find(cmp(wholeset(:),x)),subset(:),'UniformOutput',true);
        ind=arrayfun(@(x) find(cmp(wholeset(:),subset(x,:))),1:size(subset,1),'UniformOutput',true);
    catch
        %         index=cellfun(@(x) find(cmp(wholeset(:),x)),subset(:),'UniformOutput',false);
        %         nonindex=cellfun(@(x) isempty(x),index(:),'UniformOutput',true);
        %         index(nonindex)={0};
        %         index=arrayfun(@(x) find(cmp(wholeset(:),subset(x,:))),1:size(subset,1),'UniformOutput',false);
        % warning('turn to UniformOutput=false')
        % switch func2str(cmp)
        %     case {'contains','containsreverse'}
                %{
                ind=cell(size(subset,1),1);
                for i=1:size(subset,1)
                    ind{i}=find(cmp(wholeset(:),subset(i,:)));
                    if isempty(ind{i})
                        if verbose
                            warning(['missing',subset{i,:}])
                        end
                        ind{i}=0;
                    end
                    if length(ind{i})>1
                        if verbose
                            warning('match doen not unique')
                        end
                        disp(wholeset(ind{i}))
                        ind{i}=ind{i}(1);
                    end
                end
                ind=cell2mat(ind);
                %         else
                %             ind=cell(size(subset,1),1);
                %             wholeset=wholeset(:);
                %             for i=1:size(subset,1)
                %                 ind{i}=find(cmp(wholeset(i,:),subset));
                %             end
                %         end
                %}
                ind=cellfun(@(x) find(cmp(wholeset(:),x)),subset(:),'UniformOutput',false);
                ind_notfind=isemptycell(ind);
                ind(ind_notfind)={0};
                ind=cell2mat(ind);
            % case {'strcmp','strcmpi'}
            %     ind=cellfun(@(x) find(cmp(wholeset(:),x)),subset(:),'UniformOutput',false);
            %     ind_notfind=isemptycell(ind);
            %     ind(ind_notfind)={0};
            %     ind=cell2mat(ind);
            % otherwise
            %     error('not suppot this cmp')

        % end
    end
else
    ind=cellfun(@(x) find(cmp(wholeset(:),x)),subset(:),'UniformOutput',UniformOutput);
end
end






























