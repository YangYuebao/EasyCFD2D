# 采用方法A划分网格
function fieldA(x_uv::Matrix{Float64},y_uv::Matrix)
    n,m=size(x_uv)

    # 此处可以更换为更好的初场生成函数
    u_field,v_field,p_field=initialField(n,m)

    x_u,x_v=fieldDiff(x_uv,n,m)
    y_u,y_v=fieldDiff(y_uv,n,m)

    # alpha,beta,gamma的值在网格生成的时候就已经获得过了一部分，可以不用重新计算
    alpha=x_v.^2+y_v.^2
    beta=x_u.*x_v+y_u.*y_v
    gamma=x_u.^2+y_u.^2
    Ja=x_u.*y_v-x_v.*y_u

    return u_field,v_field,p_field,x_u,y_u,x_v,y_v,alpha,beta,gamma,Ja
end

# 场求导
function fieldDiff(x_uv::Matrix{Float64},n::Int64,m::Int64)
    x_u=Matrix{Float64}(undef,n,m)
    x_v=Matrix{Float64}(undef,n,m)
    x_u[:,1]=(-3*x_uv[:,1]+4*x_uv[:,2]-x_uv[:,3])/2
    x_u[:,2:end-1]=(x_uv[:,3:end]-x_uv[:,1:end-2])/2
    x_u[:,end]=(3*x_uv[:,end]-4*x_uv[:,end-1]+x_uv[:,end-2])/2
    x_v[1,:]=-(-3*x_uv[1,:]+4*x_uv[2,:]-x_uv[3,:])/2
    x_v[2:end-1,:]=-(x_uv[3:end,:]-x_uv[1:end-2,:])/2
    x_v[end,:]=-(3*x_uv[end,:]-4*x_uv[end-1,:]+x_uv[end-2,:])/2
    return x_u,x_v
end
fieldDiff(x_uv::Matrix{Float64})=fieldDiff(x_uv,size(x_uv)...)

# 生成初场
function initialField(n,m)
    u_field=zeros(n,m)
    v_field=zeros(n,m)
    p_field=zeros(n,m)
    #=
    u_field=rand(n,m).-0.5
    v_field=rand(n,m).-0.5
    p_field=rand(n,m).-0.5
    =#
    return u_field,v_field,p_field
end
