
function initialize_xy_uv(bounds::Vector{bound},n::Int64,m::Int64)
    x_uv::Matrix{Float64} = zeros(n, m)
    y_uv::Matrix{Float64} = zeros(n, m)

    temp = bounds[1].fun.(range(bounds[1].span..., length=m))
    for i = 1:m
        x_uv[end, i] = temp[i][1]
        y_uv[end, i] = temp[i][2]
    end

    temp = bounds[3].fun.(range(bounds[3].span..., length=m))
    for i = 1:m
        x_uv[1, i] = temp[end-i+1][1]
        y_uv[1, i] = temp[end-i+1][2]
    end

    temp = bounds[2].fun.(range(bounds[2].span..., length=n))
    for i = 1:n
        x_uv[i, end] = temp[end-i+1][1]
        y_uv[i, end] = temp[end-i+1][2]
    end

    temp = bounds[4].fun.(range(bounds[4].span..., length=n))
    for i = 1:n
        x_uv[i, 1] = temp[i][1]
        y_uv[i, 1] = temp[i][2]
    end

    for j = 2:m-1
        x_uv[:, j] = collect(range(x_uv[1, j], x_uv[end, j], length=n))
    end
    for j = 2:m-1
        y_uv[:, j] = collect(range(y_uv[1, j], y_uv[end, j], length=n))
    end
    return x_uv,y_uv
end

function index_generation(n::Int64,m::Int64)    
    row = Vector{Int64}(undef, 9 * m * n - 24 * (m + n) + 64) # 不变
    col = Vector{Int64}(undef, 9 * m * n - 24 * (m + n) + 64) # 不变
    n2=n-2
    # 赋值顺序：从左往右按列从上到下赋值
    # -------------------------------------------------
    # 先输入行号信息------------------------------------
    # -------------------------------------------------
    id = 1  #表示当前处理的节点在矩阵表示中的标号
    row[1:4] = ones(Int64, 4)
    id += 1
    for i = 2:n2-1
        row[6*id-7:6*id-2] = ones(Int64, 6) * id
        id += 1
    end
    row[6*id-7:6*id-4] = ones(Int64, 4) * id
    id += 1

    index_mid = 6 * n2 - 3
    for j = 2:m-3
        row[index_mid:index_mid+5] = ones(Int64, 6) * id
        id += 1
        for i = 2:n2-1
            row[(index_mid-12:index_mid-4).+9*i] = ones(Int64, 9) * id
            id += 1
        end
        row[index_mid-12+9*n2:index_mid-7+9*n2] = ones(Int64, 6) * id
        id += 1

        index_mid += (9 * n2 - 6)
    end

    row[index_mid:index_mid+3] = ones(Int64, 4) * id
    id += 1
    for i = 2:n2-1
        row[(index_mid-8:index_mid-3).+6*i] = ones(Int64, 6) * id
        id += 1
    end
    row[index_mid-8+6*n2:index_mid-5+6*n2] = ones(Int64, 4) * id

    # -------------------------------------------------
    # 输入列号信息--------------------------------------
    # -------------------------------------------------
    id = 1  #表示当前处理的节点在矩阵表示中的标号
    col[1:4] = [1, 2, 1 + n2, 2 + n2]
    id += 1
    for i = 2:n2-1
        col[6*id-7:6*id-2] = [
            id - 1,
            id,
            id + 1,
            id + n2 - 1,
            id + n2,
            id + n2 + 1
        ]
        id += 1
    end
    col[6*id-7:6*id-4] = [
        id - 1,
        id,
        id + n2 - 1,
        id + n2
    ]
    id += 1

    index_mid = 6 * n2 - 3
    for j = 2:m-3
        col[index_mid:index_mid+5] = [
            id - n2,
            id - n2 + 1,
            id,
            id + 1,
            id + n2,
            id + n2 + 1
        ]
        id += 1
        for i = 2:n2-1
            col[(index_mid-12:index_mid-4).+9*i] = [
                id - n2 - 1,
                id - n2,
                id - n2 + 1,
                id - 1,
                id,
                id + 1,
                id + n2 - 1,
                id + n2,
                id + n2 + 1
            ]
            id += 1
        end
        col[index_mid-12+9*n2:index_mid-7+9*n2] = [
            id - n2 - 1,
            id - n2,
            id - 1,
            id,
            id + n2 - 1,
            id + n2
        ]
        id += 1

        index_mid += (9 * n2 - 6)
    end

    col[index_mid:index_mid+3] = [
        id - n2,
        id - n2 + 1,
        id,
        id + 1
    ]
    id += 1
    for i = 2:n2-1
        col[(index_mid-8:index_mid-3).+6*i] = [
            id - n2 - 1,
            id - n2,
            id - n2 + 1,
            id - 1,
            id,
            id + 1,
        ]
        id += 1
    end
    col[index_mid-8+6*n2:index_mid-5+6*n2] = [
        id - n2 - 1,
        id - n2,
        id - 1,
        id
    ]
    return row,col
end

function controler(k,alpha,a,b;per=0.0)
    temp=map(x->a*(1-((x-1)/(k-1))^alpha)+b*((x-1)/(k-1))^alpha,1:k)
    #floor(Int64,per*k)
    #temp[1:floor(Int64,per*k)]=zeros
    return temp
end

function renew_coff!(n::Int64,m::Int64,alpha::Vector{Float64},beta::Vector{Float64},gamma::Vector{Float64},phi::Vector{Float64},psi::Vector{Float64},x_uv::Matrix{Float64},y_uv::Matrix{Float64})
    du=1/(m-1)
    dv=1/(n-1)
    k=(m-1)/(n-1)
    id = 1
    for j = 2:m-1
        for i = 2:n-1
            alpha[id] = ((x_uv[i+1, j] - x_uv[i-1, j]) / 2 / dv)^2 + ((y_uv[i+1, j] - y_uv[i-1, j]) / 2 / dv)^2
            beta[id] = ((x_uv[i, j+1] - x_uv[i, j-1]) * (x_uv[i+1, j] - x_uv[i-1, j]) + (y_uv[i, j+1] - y_uv[i, j-1]) * (y_uv[i+1, j] - y_uv[i-1, j])) / 4 / du / dv
            gamma[id] = ((x_uv[i, j+1] - x_uv[i, j-1]) / 2 / du)^2 + ((y_uv[i, j+1] - y_uv[i, j-1]) / 2 / du)^2            
            id += 1
        end
    end
    for j in [2,m-1]
        for i=2:n-1
            id=(i-1)+(j-2)*(n-2)
            psi[id] = -((y_uv[i+1, j] - y_uv[i-1, j]) * (y_uv[i+1, j] - 2 * y_uv[i, j] + y_uv[i-1, j]) / dv / dv + (x_uv[i+1, j] - x_uv[i-1, j]) * (x_uv[i+1, j] - 2 * x_uv[i, j] + x_uv[i-1, j]) / dv / dv) / (alpha[id]) / 2
        end
    end
    for i=1:n-2
        psi[i:(n-2):i+(m-3)*(n-2)]=controler(m-2,1.0,psi[i],psi[i+(m-3)*(n-2)])
        #psi[i:(n-2):i+(m-3)*(n-2)]=zeros(m-2)
    end
    for j = 2:m-1
        for i=[2,n-1]
            id=(i-1)+(j-2)*(n-2)
            phi[id] = -((y_uv[i, j+1] - y_uv[i, j-1]) * (y_uv[i, j+1] - 2 * y_uv[i, j] + y_uv[i, j-1]) / du / du + (x_uv[i, j+1] - x_uv[i, j-1]) * (x_uv[i, j+1] - 2 * x_uv[i, j] + x_uv[i, j-1]) / du / du) / (gamma[id]) * k / 2
        end
    end
    for i=1:m-2
        phi[(i-1)*(n-2)+1:i*(n-2)]=controler(n-2,1.0,phi[(i-1)*(n-2)+1],phi[i*(n-2)])
        #phi[(i-1)*(n-2)+1:i*(n-2)]=zeros(n-2)
    end
    return nothing
end

function get_val_bx_by!(p::Float64,q::Float64,n::Int64,m::Int64,val::Vector{Float64},bx::Vector{Float64},by::Vector{Float64},x0::Vector{Float64},y0::Vector{Float64},alpha::Vector{Float64},beta::Vector{Float64},gamma::Vector{Float64},phi::Vector{Float64},psi::Vector{Float64},x_uv::Matrix{Float64},y_uv::Matrix{Float64})
    id = 1  #表示当前处理的节点在矩阵表示中的标号
    n2=n-2
    k=(m-1)/(n-1)
    val[1:4] = [
        -2 * (alpha[id] * k * k + gamma[id]) * p,
        gamma[id]*(1 + psi[id] / 2) * q,
        (alpha[id]*(k + phi[id] / 2) * k) * q,
        (-beta[id] * k / 2) * q
    ]
    bx[id] = -(
        -beta[id] * k / 2 * x_uv[1, 1]
        + alpha[id]*(k - phi[id] / 2) * k * x_uv[2, 1]
        + beta[id] * k / 2 * x_uv[3, 1]
        + gamma[id]*(1 - psi[id] / 2) * x_uv[1, 2]
        -2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
        + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
        + beta[id] * k / 2 * x_uv[1, 3]
        + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * x0[id+n2]
        + (-beta[id] * k / 2) * (1 - q) * x0[id+n2+1]
    )
    by[id] = -(
        -beta[id] * k / 2 * y_uv[1, 1]
        + alpha[id]*( k - phi[id] / 2) * k * y_uv[2, 1]
        + beta[id] * k / 2 * y_uv[3, 1]
        + gamma[id]*(1 - psi[id] / 2) * y_uv[1, 2]
        -2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
        + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
        + beta[id] * k / 2 * y_uv[1, 3]
        + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * y0[id+n2]
        + (-beta[id] * k / 2) * (1 - q) * y0[id+n2+1]
    )
    id += 1

    #左中间
    for i = 2:n2-1
        val[6*id-7:6*id-2] = [
            gamma[id]*(1 - psi[id] / 2) * q,
            -2 * (alpha[id] * k * k + gamma[id]) * p,
            gamma[id]*(1 + psi[id] / 2) * q,
            (beta[id] * k / 2) * q,
            (alpha[id]*( k + phi[id] / 2) * k) * q,
            (-beta[id] * k / 2) * q
        ]
        bx[id] = -(
            -beta[id] * k / 2 * x_uv[id, 1]
            + alpha[id]*( k - phi[id] / 2) * k * x_uv[id+1, 1]
            + beta[id] * k / 2 * x_uv[id+2, 1]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
            -2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
            + (beta[id] * k / 2) * (1 - q) * x0[id+n2-1]
            + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * x0[id+n2]
            -beta[id] * k / 2 * (1 - q) * x0[id+n2+1]
        )
        by[id] = -(
            -beta[id] * k / 2 * y_uv[id, 1]
            + alpha[id]*( k - phi[id] / 2) * k * y_uv[id+1, 1]
            + beta[id] * k / 2 * y_uv[id+2, 1]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
            -2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
            + (beta[id] * k / 2) * (1 - q) * y0[id+n2-1]
            + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * y0[id+n2]
            -beta[id] * k / 2 * (1 - q) * y0[id+n2+1]
        )
        id += 1
    end

    #左下角
    val[6*id-7:6*id-4] = [
        gamma[id]*(1 - psi[id] / 2) * q,
        -2 * (alpha[id] * k * k + gamma[id]) * p,
        (beta[id] * k / 2) * q,
        (alpha[id]*( k + phi[id] / 2) * k) * q,
    ]
    bx[id] = -(
        -beta[id] * k / 2 * x_uv[n-2, 1]
        + alpha[id]*( k - phi[id] / 2) * k * x_uv[n-1, 1]
        + beta[id] * k / 2 * x_uv[n, 1]
        + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
        - 2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
        + gamma[id]*(1 + psi[id] / 2) * x_uv[n, 2]
        + (beta[id] * k / 2) * (1 - q) * x0[id+n2-1]
        + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * x0[id+n2]
        -beta[id] * k / 2 * x_uv[n, 3]
    )
    by[id] = -(
        -beta[id] * k / 2 * y_uv[n-2, 1]
        + alpha[id]*( k - phi[id] / 2) * k * y_uv[n-1, 1]
        + beta[id] * k / 2 * y_uv[n, 1]
        + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
        - 2 * (alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
        + gamma[id]*(1 + psi[id] / 2) * y_uv[n, 2]
        + (beta[id] * k / 2) * (1 - q) * y0[id+n2-1]
        + (alpha[id]*( k + phi[id] / 2) * k) * (1 - q) * y0[id+n2]
        -beta[id] * k / 2 * y_uv[n, 3]
    )
    id += 1

    #中间列
    index_mid = 6 * n2 - 3
    for j = 2:m-3
        # 中间的顶上
        val[index_mid:index_mid+5] = [
            alpha[id]*( k - phi[id] / 2) * k * q,
            beta[id] * k / 2 * q,
            -2(alpha[id] * k * k + gamma[id]) * p,
            gamma[id]*(1 + psi[id] / 2) * q,
            alpha[id]*( k + phi[id] / 2) * k * q,
            -beta[id] * k / 2 * q
        ]
        bx[id] = -(
            -beta[id] * k / 2 * x_uv[1, j]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
            + beta[id] * k / 2 * (1 - q) * x0[id-n2+1]
            + gamma[id]*(1 - psi[id] / 2) * x_uv[1, j+1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
            + beta[id] * k / 2 * x_uv[1, j+2]
            + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * x0[id+n2]
            -beta[id] * k / 2 * (1 - q) * x0[id+n2+1]
        )
        by[id] = -(
            -beta[id] * k / 2 * y_uv[1, j]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
            + beta[id] * k / 2 * (1 - q) * y0[id-n2+1]
            + gamma[id]*(1 - psi[id] / 2) * y_uv[1, j+1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
            + beta[id] * k / 2 * y_uv[1, j+2]
            + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * y0[id+n2]
            -beta[id] * k / 2 * (1 - q) * y0[id+n2+1]
        )
        id += 1
        # 中间列中间
        for i = 2:n2-1
            val[(index_mid-12:index_mid-4).+9*i] = [
                -beta[id] * k / 2 * q,
                alpha[id]*( k - phi[id] / 2) * k * q,
                beta[id] * k / 2 * q,
                gamma[id]*(1 - psi[id] / 2) * q,
                -2(alpha[id] * k * k + gamma[id]) * p,
                gamma[id]*(1 + psi[id] / 2) * q,
                +beta[id] * k / 2 * q,
                alpha[id]*( k + phi[id] / 2) * k * q,
                -beta[id] * k / 2 * q
            ]
            bx[id] = -(
                -beta[id] * k / 2 * (1 - q) * x0[id-n2-1]
                + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
                + beta[id] * k / 2 * (1 - q) * x0[id-n2+1]
                + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
                -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
                + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
                + beta[id] * k / 2 * (1 - q) * x0[id+n2-1]
                + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * x0[id+n2]
                - beta[id] * k / 2 * (1 - q) * x0[id+n2+1]
            )
            by[id] = -(
                -beta[id] * k / 2 * (1 - q) * y0[id-n2-1]
                + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
                + beta[id] * k / 2 * (1 - q) * y0[id-n2+1]
                + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
                -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
                + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
                + beta[id] * k / 2 * (1 - q) * y0[id+n2-1]
                + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * y0[id+n2]
                - beta[id] * k / 2 * (1 - q) * y0[id+n2+1]
            )
            id += 1
        end
        # 中间列底部
        val[index_mid-12+9*n2:index_mid-7+9*n2] = [
            -beta[id] * k / 2 * q,
            alpha[id]*( k - phi[id] / 2) * k * q,
            gamma[id]*(1 - psi[id] / 2) * q,
            -2(alpha[id] * k * k + gamma[id]) * p,
            +beta[id] * k / 2 * q,
            alpha[id]*( k + phi[id] / 2) * k * q,
        ]
        bx[id] = -(
            -beta[id] * k / 2 * (1 - q) * x0[id-n2-1]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
            + beta[id] * k / 2 * x_uv[n, j]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
            + gamma[id]*(1 + psi[id] / 2) * x_uv[n, j+1]
            + beta[id] * k / 2 * (1 - q) * x0[id+n2-1]
            + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * x0[id+n2]
            -beta[id] * k / 2 * x_uv[n, j+2]
        )
        by[id] = -(
            -beta[id] * k / 2 * (1 - q) * y0[id-n2-1]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
            + beta[id] * k / 2 * y_uv[n, j]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
            + gamma[id]*(1 + psi[id] / 2) * y_uv[n, j+1]
            + beta[id] * k / 2 * (1 - q) * y0[id+n2-1]
            + alpha[id]*( k + phi[id] / 2) * k * (1 - q) * y0[id+n2]
            -beta[id] * k / 2 * y_uv[n, j+2]
        )
        id += 1

        index_mid += (9 * n2 - 6)
    end

    # 右列右上
    val[index_mid:index_mid+3] = [
        alpha[id]*( k - phi[id] / 2) * k * q,
        beta[id] * k / 2 * q,
        -2(alpha[id] * k * k + gamma[id]) * p,
        gamma[id]*(1 + psi[id] / 2) * q,
    ]
    bx[id] = -(
        -beta[id] * k / 2 * x_uv[1, m-2]
        + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
        + beta[id] * k / 2 * (1 - q) * x0[id-n2+1]
        + gamma[id]*(1 - psi[id] / 2) * x_uv[1, m-1]
        -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
        + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
        + beta[id] * k / 2 * x_uv[1, m]
        + alpha[id]*( k + phi[id] / 2) * k * x_uv[2, m]
        -beta[id] * k / 2 * x_uv[3, m]
    )
    by[id] = -(
        -beta[id] * k / 2 * y_uv[1, m-2]
        + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
        + beta[id] * k / 2 * (1 - q) * y0[id-n2+1]
        + gamma[id]*(1 - psi[id] / 2) * y_uv[1, m-1]
        -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
        + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
        + beta[id] * k / 2 * y_uv[1, m]
        + alpha[id]*( k + phi[id] / 2) * k * y_uv[2, m]
        -beta[id] * k / 2 * y_uv[3, m]
    )
    id += 1

    # 右列中间
    for i = 2:n2-1
        val[(index_mid-8:index_mid-3).+6*i] = [
            -beta[id] * k / 2 * q,
            alpha[id]*( k - phi[id] / 2) * k * q,
            beta[id] * k / 2 * q,
            gamma[id]*(1 - psi[id] / 2) * q,
            -2(alpha[id] * k * k + gamma[id]) * p,
            gamma[id]*(1 + psi[id] / 2) * q,
        ]
        bx[id] = -(
            -beta[id] * k / 2 * (1 - q) * x0[id-n2-1]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
            + beta[id] * k / 2 * (1 - q) * x0[id-n2+1]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * x0[id+1]
            + beta[id] * k / 2 * x_uv[id-n2*(m-3), m]
            + alpha[id]*( k + phi[id] / 2) * k * x_uv[id-n2*(m-3)+1, m]
            -beta[id] * k / 2 * x_uv[id-n2*(m-3)+2, m]
        )
        by[id] = -(
            -beta[id] * k / 2 * (1 - q) * y0[id-n2-1]
            + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
            + beta[id] * k / 2 * (1 - q) * y0[id-n2+1]
            + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
            -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
            + gamma[id]*(1 + psi[id] / 2) * (1 - q) * y0[id+1]
            + beta[id] * k / 2 * y_uv[id-n2*(m-3), m]
            + alpha[id]*( k + phi[id] / 2) * k * y_uv[id-n2*(m-3)+1, m]
            -beta[id] * k / 2 * y_uv[id-n2*(m-3)+2, m]
        )
        id += 1
    end
    # 右列右下
    val[index_mid-8+6*n2:index_mid-5+6*n2] = [
        -beta[id] * k / 2 * q,
        alpha[id]*( k - phi[id] / 2) * k * q,
        gamma[id]*(1 - psi[id] / 2) * q,
        -2(alpha[id] * k * k + gamma[id]) * p,
    ]
    bx[id] = -(
        -beta[id] * k / 2 * (1 - q) * x0[id-n2-1]
        + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * x0[id-n2]
        + beta[id] * k / 2 * x_uv[n, m-2]
        + gamma[id]*(1 - psi[id] / 2) * (1 - q) * x0[id-1]
        -2(alpha[id] * k * k + gamma[id]) * (1 - p) * x0[id]
        + gamma[id]*(1 + psi[id] / 2) * x_uv[n, m-1]
        + beta[id] * k / 2 * x_uv[id-n2*(m-3), m]
        + alpha[id]*( k + phi[id] / 2) * k * x_uv[id-n2*(m-3)+1, m]
        -beta[id] * k / 2 * x_uv[id-n2*(m-3)+2, m]
    )
    by[id] = -(
        -beta[id] * k / 2 * (1 - q) * y0[id-n2-1]
        + alpha[id]*( k - phi[id] / 2) * k * (1 - q) * y0[id-n2]
        + beta[id] * k / 2 * y_uv[n, m-2]
        + gamma[id]*(1 - psi[id] / 2) * (1 - q) * y0[id-1]
        -2(alpha[id] * k * k + gamma[id]) * (1 - p) * y0[id]
        + gamma[id]*(1 + psi[id] / 2) * y_uv[n, m-1]
        + beta[id] * k / 2 * y_uv[id-n2*(m-3), m]
        + alpha[id]*( k + phi[id] / 2) * k * y_uv[id-n2*(m-3)+1, m]
        -beta[id] * k / 2 * y_uv[id-n2*(m-3)+2, m]
    )
    return val
end

