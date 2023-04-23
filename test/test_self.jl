using EasyCFD2D
begin
    k = 0
    m = 5
    n = 5
end

for j = 1:m
    for i = 1:n
        global k
        if (i == 2 || j == 2 || i == n - 1 || j == m - 1) && (i != 1) && (i != n) && (j != 1) && (j != m)
            #if (i == 2 || j == 2 || i == n - 1 || j == m - 1) 
            to_val = EasyCFD2D.to_val_index(n, m, i, j, EasyCFD2D.SecondOrderUpwind())
            k = k + sum(to_val)
            println("i=",i," j=",j," sum=",sum(to_val))
        end
    end
end

k - 24 * (m + n - 8) - 44

a= EasyCFD2D.to_val_index(n, m, 4, 3, EasyCFD2D.SecondOrderUpwind())

sum(a)