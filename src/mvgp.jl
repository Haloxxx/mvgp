module mvgp

"Finite central difference"
fdcx(f, x, yc; h=1e-5) = (f(x+h/2, yc) - f(x-h/2, yc))/h

fdcy(f, xc, y; h=1e-5) = (f(xc, y+h/2) - f(xc, y-h/2))/h

fdcx3(f, x, yc, zc; h=1e-5) = (f(x+h/2, yc, zc) - f(x-h/2, yc, zc))/h

fdcy3(f, xc, y, zc; h=1e-5) = (f(xc, y+h/2, zc) - f(xc, y-h/2, zc))/h

fdcz3(f, xc, yc, z; h=1e-5) = (f(xc, yc, z+h/2) - f(xc, yc, z-h/2))/h

function gradientDesc2(f, variables; e=1.0, β=0.01, ϵ=1e-5)
    #println("y0: ", y0)
    

    x0 = variables[1]
    y0 = variables[2]


    F0 = f(x0, y0)
    g0 = fdcx(f, x0, y0), fdcy(f, x0, y0)
    dir = g0[1]*(-1), g0[2]*(-1)
    xi = x0
    yi = y0
    #println("g0: ",g0)
    var = 1
    i = 1

    while true

        if(var == 1)
            dir = g0[1]*(-1), g0[2]*(-1)
        end
        #println("dir2: ",dir[2])
        xi_plus1 = xi + e*dir[1]
        yi_plus1 = yi + e*dir[2]

        F = f(xi_plus1, yi_plus1)
        g = fdcx(f, xi_plus1, yi_plus1), fdcy(f, xi_plus1, yi_plus1)
        #println("g: ", g)

        if F < F0 && abs(F-F0)>1e-3
            var = 1
            g0 = g
            #println("wariant1")
            xi = xi_plus1
            yi = yi_plus1
        else
            var = 2
            trans = abs(g0[1]*g0[1]+g0[2]*g0[2])
            if trans < ϵ
                break
                #println("koniec")
            end

            xi = xi_plus1 - e*dir[1]
            yi = yi_plus1 - e*dir[2]

            if(e - β > 0)
                e -= β
            end
            #println("wariant2")
        end

        if i > 1000
            break
        end
        i+=1
    end
    #println("xi: ",xi," yi: ",yi)
    #println("")
    return xi, yi
    
end

function gradientDesc3(f, variables; e=1.0, β=0.01, ϵ=1e-5)
    #println("y0: ", y0)
    

    x0 = variables[1]
    y0 = variables[2]
    z0 = variables[3]


    F0 = f(x0, y0, z0)
    g0 = fdcx3(f, x0, y0, z0), fdcy3(f, x0, y0, z0), fdcz3(f, x0, y0, z0)
    dir = g0[1]*(-1), g0[2]*(-1), g0[3]*(-1)
    xi = x0
    yi = y0
    zi = z0
    #println("g0: ",g0)
    var = 1
    i = 1

    while true

        if(var == 1)
            dir = g0[1]*(-1), g0[2]*(-1), g0[3]*(-1)
        end
        #println("dir2: ",dir[2])
        xi_plus1 = xi + e*dir[1]
        yi_plus1 = yi + e*dir[2]
        zi_plus1 = zi + e*dir[3]

        F = f(xi_plus1, yi_plus1, zi_plus1)
        g = fdcx3(f, xi_plus1, yi_plus1, zi_plus1), fdcy3(f, xi_plus1, yi_plus1, zi_plus1), fdcz3(f, x0, y0, z0)
        #println("g: ", g)

        if F < F0 && abs(F-F0)>1e-3
            var = 1
            g0 = g
            #println("wariant1")
            xi = xi_plus1
            yi = yi_plus1
            zi = zi_plus1
        else
            var = 2
            trans = abs(g0[1]*g0[1]+g0[2]*g0[2]+g0[3]+g0[3])
            if trans < ϵ
                break
                #println("koniec")
            end

            xi = xi_plus1 - e*dir[1]
            yi = yi_plus1 - e*dir[2]
            zi = zi_plus1 - e*dir[3]

            if(e - β > 0)
                e -= β
            end
            #println("wariant2")
        end

        if i > 1000
            break
        end
        i+=1
    end
    #println("xi: ",xi," yi: ",yi)
    #println("")
    return xi, yi, zi
    
end

function gradientDesc(f, variables, n; e=1.0, β=0.01, ϵ=1e-5)
    if β>=1 || β<=0
        throw("0 < β < 1 has to be fulfilled")
    end
    result = 0
    if n == 2
        result = gradientDesc2(f, variables; e = e, β = β, ϵ = ϵ)
    end

    if n == 3
        result = gradientDesc3(f, variables; e = e, β = β, ϵ = ϵ)
    end

    return result
end

export gradientDesc

end # module