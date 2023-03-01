# 计算域边界条件

m=30    #u的格数
n=30    #v的格数s

du=1/(m-1)
dv=1/(n-1)
k=(m-1)/(n-1)#dv/du

x_uv=zeros(n,m)
y_uv=zeros(n,m)

temp=bounds[1].fun.(range(bounds[1].span...,length=m))
for i=1:m
    x_uv[end,i]=temp[i][1]
    y_uv[end,i]=temp[i][2]
end

temp=bounds[3].fun.(range(bounds[3].span...,length=m))
for i=1:m
    x_uv[1,i]=temp[end-i+1][1]
    y_uv[1,i]=temp[end-i+1][2]
end

temp=bounds[2].fun.(range(bounds[2].span...,length=n))
for i=1:n
    x_uv[i,end]=temp[end-i+1][1]
    y_uv[i,end]=temp[end-i+1][2]
end

temp=bounds[4].fun.(range(bounds[4].span...,length=n))
for i=1:n
    x_uv[i,1]=temp[i][1]
    y_uv[i,1]=temp[i][2]
end

x_uv[2:n-1,2:m-1].+=sum(x_uv/(2*(m+n-2)))
x_uv[2:n-1,2:m-1]+=1e-3*rand(n-2,m-2)
y_uv[2:n-1,2:m-1].+=sum(y_uv/(2*(m+n-2)))
y_uv[2:n-1,2:m-1]+=1e-3*rand(n-2,m-2)

temp=zeros(m,n).+Inf

maxcount=20000
ep=1e-8
count=0
thisep=Inf
while thisep>ep && count<maxcount
    global count,x_uv,y_uv,temp,thisep
    temp=x_uv[:,:]
    for i=2:n-1
        for j=2:m-1
            alpha=((x_uv[i,j+1]-x_uv[i,j-1])*k/2)^2+((y_uv[i,j+1]-y_uv[i,j-1])*k/2)^2
            gamma=((x_uv[i+1,j]-x_uv[i-1,j])/2)^2+((y_uv[i+1,j]-y_uv[i-1,j])/2)^2
            beta=((x_uv[i,j+1]-x_uv[i,j-1])*(x_uv[i+1,j]-x_uv[i-1,j])+(y_uv[i,j+1]-y_uv[i,j-1])*(y_uv[i+1,j]-y_uv[i-1,j]))*k/4
            x_uv[i,j]=(alpha*(x_uv[i,j+1]+x_uv[i,j-1])*k*k+gamma*(x_uv[i+1,j]+x_uv[i-1,j])-beta*(x_uv[i+1,j+1]-x_uv[i+1,j-1]-x_uv[i-1,j+1]+x_uv[i-1,j-1])/2*k)/(alpha*k*k+gamma)/2
            y_uv[i,j]=(alpha*(y_uv[i,j+1]+y_uv[i,j-1])*k*k+gamma*(y_uv[i+1,j]+y_uv[i-1,j])-beta*(y_uv[i+1,j+1]-y_uv[i+1,j-1]-y_uv[i-1,j+1]+y_uv[i-1,j-1])/2*k)/(alpha*k*k+gamma)/2
        end
    end
    count+=1
    thisep=maximum(abs.(temp-x_uv))
    println(count," log10(ep)=",log10(thisep))
end

myplot(x_uv,y_uv,m,n)

if count>=maxcount
    @warn "Not convergence!"
end

function myplot(x_uv,y_uv,m,n)
    plt = plot(
        legend=:none,
        xlim=(0, 3),
        ylim=(0, 2)
    )
    for i = 1:n
        plot!(plt, x_uv[i, :][:], y_uv[i, :][:], linecolor=:blue)
    end

    for j = 1:m
        plot!(plt, x_uv[:, j][:], y_uv[:, j][:], linecolor=:blue)
    end
    return plt
end
