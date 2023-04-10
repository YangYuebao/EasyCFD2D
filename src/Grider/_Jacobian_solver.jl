using SparseArrays
using EquationsSolver
# 计算域边界条件
begin
    m = 35    #u的格数
    n = 40    #v的格数
    m2 = m - 2 
    n2 = n - 2
    du = 1 / (m - 1)
    dv = 1 / (n - 1)
    k = (m - 1) / (n - 1)#dv/du
    dims=floor(Int64,(m*n)^(0.3))
end

x_uv,y_uv=initialize_xy_uv(bounds,n,m)

row,col=index_generation(n,m)

begin # new variables
    val = Vector{Float64}(undef, 9 * m * n - 24 * (m + n) + 64) # 系数矩阵
    alpha = Vector{Float64}(undef, (m - 2) * (n - 2))       # 
    beta = Vector{Float64}(undef, (m - 2) * (n - 2))        # 
    gamma = Vector{Float64}(undef, (m - 2) * (n - 2))       # 
    phi = Vector{Float64}(undef, (m - 2) * (n - 2))         # 
    psi = Vector{Float64}(undef, (m - 2) * (n - 2))         # 
    bx = Vector{Float64}(undef, (m - 2) * (n - 2))          # x非齐次项
    by = Vector{Float64}(undef, (m - 2) * (n - 2))          # y非齐次项
    x0 = reshape(x_uv[2:n-1, 2:m-1], :)                       # x估计值
    y0 = reshape(y_uv[2:n-1, 2:m-1], :)                       # y估计值
    p = 1.0   # 主对角线松弛系数
    q = 1.0   # 相邻节点松弛系数s
    ep=1
    maxep=1e-8
    mycount=0
    maxmycount=100000
end
# -------------------------------------------------
# 获取各种系数--------------------------------------
# -------------------------------------------------

#begin 
    #for myasdf=1:100
    while ep>maxep && ep<5 && mycount<maxmycount    
    global ep,mycount
    renew_coff!(n,m,alpha,beta,gamma,phi,psi,x_uv,y_uv)

    # -------------------------------------------------
    # 输入系数信息--------------------------------------
    # -------------------------------------------------
    #左上角
    A=get_A_bx_by!(n,m,val,bx,by,alpha,beta,gamma,phi,psi,x_uv,y_uv)
    
    tempx = x0[:]
    tempy = y0[:]

    #=
    lpx = LinearProblem(A, bx, x0)
    lpy = LinearProblem(A, by, y0)
    x0 = solve(lpx, GMRESM(); m=dims, maxiter=10+mycount,abstol=maxep,warning=true)
    y0 = solve(lpy, GMRESM(); m=dims, maxiter=10+mycount,abstol=maxep,warning=true)
    =#

    x0=A\bx
    y0=A\by
    
    relax=0.1
    x0=relax*x0+(1-relax)*tempx
    y0=relax*y0+(1-relax)*tempy
    

    ep = maximum(abs.(x0 - tempx)) + maximum(abs.(y0 - tempy))

    x_uv[2:n-1, 2:m-1] = reshape(x0, n2, m2)[:,:]
    y_uv[2:n-1, 2:m-1] = reshape(y0, n2, m2)[:,:]
    
    if mycount%10==0
        println(mycount," log10(ep)=",log10(ep))
        display(myplot(x_uv, y_uv, m, n))
    end
    mycount+=1
end


if mycount >= maxmycount
    @warn "Not convergence!"
end
plt=myplot(x_uv, y_uv, m, n)
savefig(plt,"pic.png")