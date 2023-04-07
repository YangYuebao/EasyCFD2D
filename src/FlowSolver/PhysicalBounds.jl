#=
function fieldProblemB(x_uv::Matrix{Float64},y_uv::Matrix)
    x_field=(x_uv[1:end-1,1:end-1]+x_uv[1:end-1,2:end]+x_uv[2:end,1:end-1]+x_uv[2:end,2:end])/4
    y_field=(y_uv[1:end-1,1:end-1]+y_uv[1:end-1,2:end]+y_uv[2:end,1:end-1]+y_uv[2:end,2:end])/4
    n,m=size(x_field)

    # 此处可以更换为更好的初场生成函数
    u_field=zeros(n,m)
    v_field=zeros(n,m)
    p_field=zeros(n,m)

    x_u=(x_uv[1:end-1,2:end]-x_uv[1:end-1,1:end-1]+x_uv[2:end,2:end]-x_uv[2:end,1:end-1])/2
    y_u=(y_uv[1:end-1,2:end]-y_uv[1:end-1,1:end-1]+y_uv[2:end,2:end]-y_uv[2:end,1:end-1])/2
    x_v=(x_uv[2:end,2:end]+x_uv[2:end,1:end-1]-x_uv[1:end-1,2:end]-x_uv[1:end-1,1:end-1])/2
    y_v=(y_uv[2:end,2:end]+y_uv[2:end,1:end-1]-y_uv[1:end-1,2:end]-y_uv[1:end-1,1:end-1])/2
    alpha=x_v.^2+y_v.^2
    beta=x_u.*x_v+y_u.*y_v
    gamma=x_u.^2+y_u.^2
    Ja=x_u.*y_v-x_v.*y_u

    return x_field,y_field,u_field,v_field,p_field,x_u,y_u,x_v,y_v,alpha,beta,gamma,Ja
end
=#

function fieldA(x_uv::Matrix{Float64},y_uv::Matrix)
    n,m=size(x_field)

    # 此处可以更换为更好的初场生成函数
    u_field=zeros(n,m)
    v_field=zeros(n,m)
    p_field=zeros(n,m)

    x_u=Matrix{Float64}(undef,n,m)
    y_u=Matrix{Float64}(undef,n,m)
    x_v=Matrix{Float64}(undef,n,m)
    y_v=Matrix{Float64}(undef,n,m)

    x_u[:,1]=(-3*x_uv[:,1]+4*x_uv[:,2]-x_uv[:,3])/2
    x_u[:,2:end-1]=(x_uv[:,3:end]-x_uv[:,1:end-2])/2
    x_u[:,end]=(3*x_uv[:,end-2]-4*x_uv[:,end-1]+x_uv[:,end-2])/2

    y_u[:,1]=(-3*y_uv[:,1]+4*y_uv[:,2]-y_uv[:,3])/2
    y_u[:,2:end-1]=(y_uv[:,3:end]-y_uv[:,1:end-2])/2
    y_u[:,end]=(3*y_uv[:,end-2]-4*y_uv[:,end-1]+y_uv[:,end-2])/2

    x_v[1,:]=(-3*x_uv[1,:]+4*x_uv[2,:]-x_uv[3,:])/2
    x_v[2:end-1,:]=(x_uv[3:end,:]-x_uv[1:end-2,:])/2
    x_v[end,:]=(3*x_uv[end-2,:]-4*x_uv[end-1,:]+x_uv[end-2,:])/2

    y_v[1,:]=(-3*y_uv[1,:]+4*y_uv[2,:]-y_uv[3,:])/2
    y_v[2:end-1,:]=(y_uv[3:end,:]-y_uv[1:end-2,:])/2
    y_v[end,:]=(3*y_uv[end-2,:]-4*y_uv[end-1,:]+y_uv[end-2,:])/2

    alpha=x_v.^2+y_v.^2
    beta=x_u.*x_v+y_u.*y_v
    gamma=x_u.^2+y_u.^2
    Ja=x_u.*y_v-x_v.*y_u

    return u_field,v_field,p_field,x_u,y_u,x_v,y_v,alpha,beta,gamma,Ja
end

# 生成初场
function initialField()

end

# 生成源项
function sorceTerm()
    
end