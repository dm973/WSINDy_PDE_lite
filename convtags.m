function tags = convtags(varargin)
   defaultExp = [];
   defaultMon = [];
   defaultLag = [];
   defaultSing = [];
   defaultGauss = [];
   defaultSingeps = 10^-2;
   defaultsvdtol = [];
   defaultutagin = [];
   defaultutagout = [];
   defaultpsitags = [];

   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && all(x >= 0);
   addParameter(p,'Exp',defaultExp,validScalarPosNum);
   addParameter(p,'Mon',defaultMon,validScalarPosNum);
   addParameter(p,'Lag',defaultLag,validScalarPosNum);
   addParameter(p,'Sing',defaultSing,@(x)validScalarPosNum(-x));
   addParameter(p,'Gauss',defaultGauss);
   addParameter(p,'Singeps',defaultSingeps);
   addParameter(p,'svdtol',defaultsvdtol);
   addParameter(p,'utagin',defaultutagin);
   addParameter(p,'utagout',defaultutagout);
   addParameter(p,'psitags',defaultpsitags);
   parse(p,varargin{:});
   
   utagin = p.Results.utagin;
   utagout = p.Results.utagout;
   psitags = p.Results.psitags;  
   dimx=size(psitags,2)-1;
    
   tags = {};

   for jj=1:size(utagin,1)       
    
       utagin_temp = utagin(jj,:);
       utagout_temp = utagout(jj,:);
       if any(psitags(jj,:)<0)
            psitags_temp = -psitags(jj,:);
            f = {@(h,a,r,varargin) h(r,a)};
       else
            psitags_temp = [eye(dimx)*psitags(1,jj) zeros(dimx,1)];
            f = cellfun(@(i) @(h,a,r,varargin) h(r,a)./max(r,eps).*varargin{min(i,end)},  num2cell(1:dimx), 'uni',0);
       end
       Expparams = p.Results.Exp;
       ExpL = length(Expparams);
       Expcell = cell(dimx*ExpL,1);
       Monparams = p.Results.Mon;
       MonL = length(Monparams); 
       Moncell = cell(dimx*MonL,1);
       Lagparams = p.Results.Lag;
       LagL = length(Lagparams); 
       Lagcell = cell(dimx*LagL,1);
       Singparams = p.Results.Sing;
       SingL = length(Singparams);
       Gaussparams = p.Results.Gauss;
       GaussL = size(Gaussparams,2);
       Gausscell = cell(dimx*GaussL,1);
       Singcell = cell(dimx*SingL,1);
       svdtol = p.Results.svdtol;


       exptag = 'e';
       montag = 'mon';
       singtag = 'sing';
       lagtag = 'lag';
       gausstag = 'gauss';

%        for i=1:ExpL
%            for j=1:length(f)
%                Ktag_temp = ['*',exptag,num2str(Expparams(i))];
%                ui = build_str_tags([utagin_temp psi_tags_temp(j,:)],length(utagin_temp),length(psi_tags_temp(j,:)));
%                uo = build_str_tags([utagout_temp zeros(1,length(psi_tags_temp(j,:)))],length(utagin_temp),length(psi_tags_temp(j,:)));
%                tag=strrep(ui{1},'_',[Ktag_temp,'))_']);
%                tag=['(',uo{1}(1:strfind(uo{1},'_')-1),'(',tag];
%                Expcell{length(f)*i-(length(f)-j)} = {utagin_temp, @(varargin) f{j}(h,ndnormsum(varargin{:}),varargin{:}),psitags_temp(j,:),utagout_temp,svdtol,strrep(tag,'*',['*',dimtags(j)])};
%            end
%        end

       Expcell = buildKinfo(@(r,a) a*exp(-a*r),exptag,Expparams,f,utagin_temp,psitags_temp,utagout_temp,svdtol);

       Moncell = buildKinfo(@(r,a)r.^a,montag,Monparams,f,utagin_temp,psitags_temp,utagout_temp,svdtol);

       Lagcell = buildKinfo(@(r,a) horn(LaguerrePoly(a),r),lagtag,Lagparams,f,utagin_temp,psitags_temp,utagout_temp,svdtol);

       if GaussL>0
           Gausscell = buildKinfo(@(r,a) exp(-(r-a).^2/2/Gaussparams(2,1)^2),gausstag,Gaussparams(:,1),f,utagin_temp,psitags_temp,utagout_temp,svdtol);
       end
       Singcell = buildKinfo(@(r,a) (a~=0)*max(r,p.Results.Singeps).^a + (a==0).*log(max(r,p.Results.Singeps)),singtag,Singparams,f,utagin_temp,psitags_temp,utagout_temp,svdtol);
       
%        for i=1:MonL
%            Ktag_temp = ['*',montag,num2str(Monparams(i))];
%            if dimx == 1
%                Moncell{i} = {utagin_temp, @(x,t) f(@(x)x.^Monparams(i),abs(x)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convmonx',num2str(Monparams(i))]};
%            elseif dimx == 2
%                Moncell{2*i-1} = {utagin_temp, @(x,y,t) f(@(r)r.^Monparams(i),hypot(x,y)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convmonx',num2str(Monparams(i))]};
%                Moncell{2*i} = {utagin_temp, @(x,y,t) f(@(r)r.^Monparams(i),hypot(x,y)).*y,psitags_temp(2,:),utagout_temp,svdtol,['convmony',num2str(Monparams(i))]};
%            end
%        end
%        
%        for i=1:LagL
%            Ktag_temp = ['*',lagtag,num2str(Lagparams(i))];
%            c = LaguerrePoly(Lagparams(i));
%            if dimx == 1
%                Lagcell{i} = {utagin_temp, @(x,t) f(@(x) horn(c,x),abs(x)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convlagx',num2str(Lagparams(i))]};
%            elseif dimx == 2
%                Lagcell{2*i-1} = {utagin_temp, @(x,y,t) f(@(r) horn(c,r),hypot(x,y)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convmonx',num2str(Lagparams(i))]};
%                Lagcell{2*i} = {utagin_temp, @(x,y,t) f(@(r) horn(c,r),hypot(x,y)).*y,psitags_temp(2,:),utagout_temp,svdtol,['convmony',num2str(Lagparams(i))]};
%            end
%        end
% 
%        for i=1:GaussL
%            Ktag_temp = ['*',gausstag,num2str(Gaussparams(1,i))];
%            if dimx == 1
%                Gausscell{i} = {utagin_temp, @(x,t) f(@(x)exp(-(x-Gaussparams(1,i)).^2/2/Gaussparams(2,i)^2),abs(x)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convgaussx',num2str(Gaussparams(1,i))]};
%            elseif dimx == 2
%                Gausscell{2*i-1} = {utagin_temp, @(x,y,t) f(@(r)exp(-(r-Gaussparams(1,i)).^2/2/Gaussparams(2,i)^2),hypot(x,y)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convgaussx',num2str(Gaussparams(1,i))]};
%                Gausscell{2*i} = {utagin_temp, @(x,y,t) f(@(r)exp(-(r-Gaussparams(1,i)).^2/2/Gaussparams(2,i)^2),hypot(x,y)).*y,psitags_temp(1,:),utagout_temp,svdtol,['convgaussy',num2str(Gaussparams(1,i))]};
%            end
%        end
% 
%        for i=1:SingL
%            Ktag_temp = ['*',singtag,num2str(Singparams(i))];
%            if dimx == 1
%                if Singparams(i)<0
%                    Singcell{i} = {utagin_temp, @(x,t) f(@(x)(max(x,p.Results.Singeps)).^(Singparams(i)),abs(x)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convsingx',num2str(Singparams(i))]};
%                elseif Singparams(i)==0
%                    Singcell{i} = {utagin_temp, @(x,t) f(@(x)log(max(x,p.Results.Singeps)),abs(x)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convsingx',num2str(Singparams(i))]};
%                end
%            elseif dimx == 2
%                if Singparams(i)<0
%                    Singcell{2*i-1} = {utagin_temp, @(x,y,t) f(@(r)(max(r,p.Results.Singeps)).^(Singparams(i)),hypot(x,y)).*x,psitags_temp(1,:),utagout_temp,svdtol,['convsingx',num2str(Singparams(i))]};
%                    Singcell{2*i} = {utagin_temp, @(x,y,t) f(@(r)(max(r,p.Results.Singeps)).^(Singparams(i)),hypot(x,y)).*y,psitags_temp(2,:),utagout_temp,svdtol,['convsingy',num2str(Singparams(i))]};
%                elseif Singparams(i)==0
%                    Singcell{2*i-1} = {utagin_temp, @(x,y,t) f(@(r)log(max(r,p.Results.Singeps)),hypot(x,y)).*x,psitags_temp(1,:),utagout_temp,svdtol,['*singx',num2str(Singparams(i))]};
%                    Singcell{2*i} = {utagin_temp, @(x,y,t) f(@(r)log(max(r,p.Results.Singeps)),hypot(x,y)).*y,psitags_temp(2,:),utagout_temp,svdtol,['*singy',num2str(Singparams(i))]};
%                end
%            end
%        end 
       tags = [tags;[Expcell;Moncell;Lagcell;Singcell;Gausscell]];
   end
end


function out = buildKinfo(h,htag,params,f,utagin_temp,psitags_temp,utagout_temp,svdtol)
    dimtags='xyz';
    out = {};
    for i=1:length(params)
        for j=1:length(f)
            Ktag_temp = ['conv',dimtags(j),'_',htag,num2str(params(i))];
            ui = build_str_tags([utagin_temp psitags_temp(j,:)],length(psitags_temp(j,:)),length(utagin_temp));
            uo = build_str_tags([utagout_temp zeros(1,length(psitags_temp(j,:)))],length(psitags_temp(j,:)),length(utagin_temp));
            tag=['(',uo{1}(1:strfind(uo{1},'_')-1),'(',strrep(ui{1},'_',[Ktag_temp,'))_'])];
            out{end+1} = {utagin_temp, @(varargin) f{j}(h,params(i),ndnormsum(varargin{:}),varargin{:}),psitags_temp(j,:),utagout_temp,svdtol,tag};
        end
    end
    out=out(:);
end
