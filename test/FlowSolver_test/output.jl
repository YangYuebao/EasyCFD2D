function output(A)
    file=open("test1.csv","w")
    n,m=size(A)
    for i=1:n
        for j=1:m
            write(file,string(A[i,j],","))
        end
        write(file,"\n")
    end
    close(file)
end