function    Distribution = cvxregsolver(Method,Signal,L,Kernel,RegParam)

Dimension = length(Signal);

switch Method
    case 'tikhonov'
        cvx_begin quiet
            variable Distribution(Dimension)
            minimize( 1/2*square_pos(norm(Kernel*Distribution - Signal,'fro')) + 1/2*square_pos(norm(RegParam*L*Distribution,'fro')))
            subject to
                Distribution >= 0;
        cvx_end
    case 'tv'
        cvx_begin quiet
            variable Distribution(Dimension)
            clear VectorialNorm
            Auxiliary = L*Distribution;
            % Vectorial norm has not implemented until MATLAB R2017b, do manually
            for i=1:length(Auxiliary)
                %Ignore preallocation warning to ensure 'cvx convex expression' class
                VectorialNorm(i) =  norm([Auxiliary(i) 1e-12]);
            end
            minimize( 1/2*square_pos(norm(Kernel*Distribution - Signal,'fro')) + RegParam^2*sum(VectorialNorm) )
            subject to
                Distribution >= 0;
        cvx_end
end


end