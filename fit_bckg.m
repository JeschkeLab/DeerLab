function [bckg,handback]=fit_bckg(handles,texp,t_fit,data_fit)
%
% Fits background decay according to the model selected by
% handles.background
% 0  fractal, n variable,  exp(-k*t^(n/3))
% 1  n-dimensional, n fixed, exp(-k*t^(n/3))
% 2  three-dimensional, exp(-k*t)
% 3  polynomial
% 4  user-defined function in handles.bckg_fct or numerical background in handles.bckg_data
% 5  homogeneous background 
%
% texp      full time range
% t_fit     time range for background fit
% data_fit  time-domain data over t_fit

back_model = 1;
man_bckg_flag=get(handles.manual_bckg,'Value');

polynomialBckg = get(handles.bckg_poly,'Value');
if polynomialBckg, back_model = 3; end
homogeneousBg = get(handles.bckg_homogeneous,'Value');
if homogeneousBg
    fitDimension = get(handles.bckg_fit_dim,'Value');
    if fitDimension
        back_model=0;
    else
      if handles.hom_dim~=3
        back_model=1;
      else
        back_model=2;
      end
    end
    if man_bckg_flag
        back_model=5;
    end
end
if back_model == 1
    set(handles.bckg_fit_dim,'Value',0);
    set(handles.manual_bckg,'Value',0);
end
userBckg=get(handles.bckg_exp,'Value');
if userBckg
    back_model=4;
    set(handles.density_text,'String','Rel. dens.');
else
    set(handles.density_text,'String','Density');
end    
    
poly=polyfit(t_fit,log(data_fit),1); % linear fit of logarithm
v0=[-poly(1) 1];
switch back_model
    case 0
        v0=[-poly(1) 1 3];
        v1=fminsearch(@rms_stretched,v0,[],t_fit,data_fit,handles);
        % bckg=decay_stretched(v1,texp);
        % distr=handles.Pake_r.^(v1(3)-1);
        % distr=0.01*distr/sum(distr);
        % logB=distr*handles.Pake_kernel;
        bckg=decaynD(v1(1:2),texp,v1(3));
        density=v1(1);
        fdim=v1(3);
        pstr=sprintf('%5.2f',fdim);
        set(handles.bckg_dim_edit,'String',pstr);
        handles.hom_dim=fdim;
        handles.bckg_dens=density;
    case 1
        distr=handles.Pake_r.^(handles.hom_dim-1);
        distr=0.01*distr/sum(distr);
        logB=distr*handles.Pake_kernel;
        v1=fminsearch(@rmsnD,v0,[],t_fit,data_fit,logB,handles);
        bckg=decaynD(v1,texp,handles.hom_dim);
        density=v1(1);
        handles.bckg_dens=density;
    case 2
		bckg=exp(polyval(poly,abs(texp))); % background is exponential of that
        density=v0(1);
        handles.bckg_dens=density;
    case 3
        poly=polyfit(t_fit,log(data_fit),handles.poly_order);
        bckg=exp(polyval(poly,abs(texp)));
        handles.polynomial=poly;
        density=-1;
    case 4
        bckg0=exp(polyval(handles.polynomial,abs(t_fit)));
        v0=[1 0.8];
        v1=fminsearch(@rms_ubckg,v0,[],data_fit,bckg0);
        bckg1=exp(polyval(handles.polynomial,abs(texp)));
        bckg=v1(2)*exp(v1(1)*log(bckg1));
        density=v1(1)/handles.calib_density;
        handles.bckg_dens=density;
    case 5
        %targ=texp.^(handles.hom_dim/3);
        v1(1)=handles.man_k;
        v1(2)=1;
        bckg=decaynD(v1,texp,handles.hom_dim);
        bckg=(1-handles.man_depth)*bckg;
        % bckg=(1-handles.man_depth)*exp(-dens*targ);
        density=handles.man_k;
end
% figure(13); clf;
% hold on;
% plot(t_fit,td_fit,'k');
% plot(t_fit,bckg0,'g');
% plot(texp,bckg,'r');
% figure(14);
% bckg=real(bckg);
pstr=sprintf('%6.3f',density*handles.calib_density);
if density>=0
	set(handles.bckg_density,'String',pstr);
else
	set(handles.bckg_density,'String','n.a.');
end
handles.density_value=density*handles.calib_density;

handback=handles;

bckg=real(bckg);
