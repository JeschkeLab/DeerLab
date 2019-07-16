function    cvx_distribution = cvxregsolver(Method,Signal,L,Kernel,RegParam,NonNegConst)

Dimension = length(Signal);

switch Method
    case 'tikhonov'
        cvx_begin quiet
            variable cvx_distribution(Dimension)
            minimize( 1/2*square_pos(norm(Kernel*cvx_distribution - Signal,'fro')) + 1/2*square_pos(norm(RegParam*L*cvx_distribution,'fro')))
            if NonNegConst
            subject to
                cvx_distribution >= 0;
            end
        cvx_end
    case 'tv'
        cvx_begin quiet
            variable cvx_distribution(Dimension)
            clear VectorialNorm
            Auxiliary = L*cvx_distribution;
            % Vectorial norm has not implemented until MATLAB R2017b, do manually
            for i=1:length(Auxiliary)
                %Ignore preallocation warning to ensure 'cvx convex expression' class
                VectorialNorm(i) =  norm([Auxiliary(i) 1e-12]);
            end
            minimize( 1/2*square_pos(norm(Kernel*cvx_distribution - Signal,'fro')) + RegParam^2*sum(VectorialNorm) )
            if NonNegConst
                subject to
                cvx_distribution >= 0;
            end
        cvx_end
end


end