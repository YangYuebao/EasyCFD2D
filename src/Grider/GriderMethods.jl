# 高斯赛德尔迭代拉普拉斯法
function gs!(m::Int64,n::Int64,x_uv::Matrix{Float64}, y_uv::Matrix{Float64},k::Float64; maxcount::Int64=10000, ep::Float64=1e-8)
    temp = zeros(m, n) .+ Inf
    mycount = 0
    thisep = Inf
    while thisep > ep && mycount < maxcount
        temp = x_uv[:, :]
        for i = 2:n-1
            for j = 2:m-1
                alpha = ((x_uv[i, j+1] - x_uv[i, j-1]) * k / 2)^2 + ((y_uv[i, j+1] - y_uv[i, j-1]) * k / 2)^2
                gamma = ((x_uv[i+1, j] - x_uv[i-1, j]) / 2)^2 + ((y_uv[i+1, j] - y_uv[i-1, j]) / 2)^2
                beta = ((x_uv[i, j+1] - x_uv[i, j-1]) * (x_uv[i+1, j] - x_uv[i-1, j]) + (y_uv[i, j+1] - y_uv[i, j-1]) * (y_uv[i+1, j] - y_uv[i-1, j])) * k / 4
                x_uv[i, j] = (gamma * (x_uv[i, j+1] + x_uv[i, j-1]) * k * k + alpha * (x_uv[i+1, j] + x_uv[i-1, j]) - beta * (x_uv[i+1, j+1] - x_uv[i+1, j-1] - x_uv[i-1, j+1] + x_uv[i-1, j-1]) / 2 * k) / (alpha + gamma * k * k) / 2
                y_uv[i, j] = (gamma * (y_uv[i, j+1] + y_uv[i, j-1]) * k * k + alpha * (y_uv[i+1, j] + y_uv[i-1, j]) - beta * (y_uv[i+1, j+1] - y_uv[i+1, j-1] - y_uv[i-1, j+1] + y_uv[i-1, j-1]) / 2 * k) / (alpha + gamma * k * k) / 2
            end
        end
        mycount += 1
        thisep = maximum(abs.(temp - x_uv))
        println(mycount, " log10(ep)=", log10(thisep))
    end
    return x_uv,y_uv
end

function GSGrider(m::Int64,n::Int64,bounds::Vector{bound};maxcount::Int64=10000,maxep::Float64=1e-8)
    k = (m - 1) / (n - 1)#dv/du
    x_uv,y_uv=initialize_xy_uv(bounds,n,m)
    gs!(m,n,x_uv, y_uv,k; maxcount=maxcount, ep=maxep)
    return x_uv,y_uv
end

# 雅克比迭代边界正交改良法
function JacobianGrider(m::Int64,n::Int64,bounds::Vector{bound};maxcount::Int64=10000,maxep::Float64=1e-3,relax::Float64=1.0,displayStep::Int64=20)
    x_uv,y_uv=initialize_xy_uv(bounds,n,m)
    row,col=index_generation(n,m)

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
    mycount=0

    while ep>maxep && ep<5 && mycount<maxcount    
        renew_coff!(n,m,alpha,beta,gamma,phi,psi,x_uv,y_uv)
        get_val_bx_by!(p,q,n,m,val,bx,by,x0,y0,alpha,beta,gamma,phi,psi,x_uv,y_uv)
        A = sparse(row, col, val)

        tempx = x0[:]
        tempy = y0[:]
    
        x0=A\bx
        y0=A\by
        x0=relax*x0+(1-relax)*tempx
        y0=relax*y0+(1-relax)*tempy        
    
        ep = maximum(abs.(x0 - tempx)) + maximum(abs.(y0 - tempy))
    
        x_uv[2:n-1, 2:m-1] = reshape(x0, n-2, m-2)[:,:]
        y_uv[2:n-1, 2:m-1] = reshape(y0, n-2, m-2)[:,:]
        
        if mycount % displayStep == 0
            println(mycount," log10(ep)=",log10(ep))
            #display(myplot(x_uv, y_uv, m, n))
        end
        mycount+=1
    end
    println(mycount-1," log10(ep)=",log10(ep))
    return x_uv,y_uv
end