module Ising

export Estado
export energias_t
export magnetizaciones_t
export energiapromedio_beta
export magnetizacionpromedio_beta

type Estado
    sigma::Array{Float64,2}
    E::Real
    M::Real
end

function conf_aleatoria(n::Int64,m::Int64,p=0.5)
    configuracion=ones(n,m)
    for i in 1:n
        for j in 1:m
            if rand()<=p
                configuracion[i,j]=-1
            end
        end
    end
    return configuracion
end

function flip_one(A::Array{Float64,2},i::Int64,j::Int64)
    A[i,j]*=-1
    A
end

function energia_ij(configuracion::Array{Float64,2},n::Int64,m::Int64,i::Int64,j::Int64)
    -configuracion[i,j]*(configuracion[mod1(i-1,n),j]+configuracion[mod1(i+1,n),j]+
        configuracion[i,mod1(j-1,m)]+configuracion[i,mod1(j+1,m)])/2
end

function α(configuracion::Array{Float64,2},beta::Float64,n::Int64,m::Int64,i::Int64,j::Int64)
    energia = energia_ij(configuracion,n,m,i,j)
    if energia <=0.0
        delta_energia=-2*energia
        return e^(-beta*delta_energia)
    else
        return 1.0
        #energia_ij(configuracion_new,n,m,i,j)-energia_ij(configuracion,n,m,i,j)
        #min(1,e^(-beta*delta_energia))
    end
end

function aceptar(configuracion::Array{Float64,2},beta::Float64,n::Int64,m::Int64,i::Int64,j::Int64)
    alpha=α(configuracion,beta,n,m,i,j)
    #@show alpha
    if rand()<=alpha
        return 1
    else
        return 0
    end
end

function energia_total(configuracion::Array{Float64,2},n::Int64,m::Int64)
    out=0.0
    for i in 1:n
        for j in 1:m
            out+=-configuracion[i,j]*(configuracion[mod1(i-1,n),j]+configuracion[mod1(i+1,n),j]+
            configuracion[i,mod1(j-1,m)]+configuracion[i,mod1(j+1,m)])
        end
    end
    out/2
end

function energias_t(beta,n::Int64,m::Int64,t=100)
    out=zeros(t+1)
    config=conf_aleatoria(n,m)
    out[1]=energia_total(config,n,m)
    for tiempo in 1:t
        i,j=rand(1:n),rand(1:m)
        if aceptar(config,beta,n,m,i,j)<1
        out[tiempo+1]=out[tiempo]
        else
        out[tiempo+1]=out[tiempo]-2*energia_ij(config,n,m,i,j)
        config=flip_one(config,i,j)
        end
    end
    out
end

function magnetizaciones_t(beta,n::Int64,m::Int64,t)
    out=zeros(t+1)
    config=conf_aleatoria(n,m)
    out[1]=sum(config)
    for tiempo in 1:t
        i,j=rand(1:n),rand(1:m)
        if aceptar(config,beta,n,m,i,j)<1
        out[tiempo+1]=out[tiempo]
        else
        out[tiempo+1]=out[tiempo]-2*config[i,j]
        config=flip_one(config,i,j)
        end
    end
    out
end

function energiapromedio_temp(temp_max::Float64,n::Int64,m::Int64,t::Float64,paso::Float64=0.1)
    temperaturas=[0:paso:temp_max]
    maximo=length(temperaturas)
    out=zeros(maximo)
    for i in 1:maximo
        out[i]=mean(energias_t(1/temperaturas[i],n,m,t))
    end
    out
end

function magnetizacionpromedio_temp(temp_max::Float64,n::Int64,m::Int64,t::Float64,paso::Float64=0.1)
    temperaturas=[0:paso:temp_max]
    maximo=length(temperaturas)
    out=zeros(maximo)
    for i in 1:maximo
        out[i]=mean(magnetizaciones_t(1/temperaturas[i],n,m,t))
    end
    out
end

function microestados(n :: Int64, m :: Int64)
    N=n*m
    if N>16 return 0 end
    N2=big(2)^N
    out=zeros(N2,N)
    for i in 1:N2-1
        bini=bin(i)
        out[i+1,N-length(bini)+1:N]=(int(split(bini,"")).*2)-1
    end
    out
end

function particion_T(T :: Float64, configuraciones :: Array{Float64,2}, n :: Int64, m :: Int64)
    out=0.0
    for i in length(configuraciones[:,1])
        out+=e^(-energia_total(configuraciones[i,:],n,m)/T)
    end
    out
end

end
