module Ising

export MicroEstado, edo_aleatorio, voltea_espin!, energia_total, energia_ij
export propone_cambio, paso_monte_carlo, simulacion_monte_carlo
import Base.show

type MicroEstado
    σ::Array{Int,2}
	#Vamos a suponer que todas las configuraciones son cuadradas
    L::Int
end

show(io::IO, m::MicroEstado) = print(io, m.σ)

function edo_aleatorio(L::Int)
    σ = ones(Int, (L,L))
    for i in 1:L^2
        if rand() <= 0.5
			σ[i] = -1
        end
	end
    MicroEstado(σ,L)
end

function voltea_espin!(m::MicroEstado, i::Int, j::Int)
    m.σ[i,j] *= -1
end

function energia_total(m::MicroEstado)
	out = 0.
    for i in 1:m.L, j in m.L
		out -= m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
    end
    out/2
end

function energia_ij(m::MicroEstado, i::Int, j::Int)
	- m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
end

function propone_cambio(m::MicroEstado, β::Float64)
    i, j = rand(1:m.L), rand(1:m.L)  #Es más rápido que rand(1:m.L, 2)
	ΔE = -2*energia_ij(m, i, j)

	ΔE, i, j
end

##VA A REGRESAR A CASA PASO EL ESTADO (ARAMIS) O NO (NACHO)?

function paso_monte_carlo(m::MicroEstado, β::Float64)
	ΔE, i, j = propone_cambio(m, β)

	#El parámetro 	de aceptación
	α = min(1., e^(-β*ΔE))

    if rand() < α
        voltea_espin!(m, i, j)
        return ΔE
    else
        return 0
    end
end

function simulacion_monte_carlo(L::Int, T::Float64, num_pasos::Int)
	β = 1/T
	m = edo_aleatorio(L)

	for i in 1:num_pasos
		paso_monte_carlo(m,β)
	end

	m
end









function energias_t(beta,n::Int,m::Int,steps=100)
    out=zeros(steps+1)
    config_old=conf_aleatoria(n,m)
    out[1]=energia_total(config_old,n,m)
    for tiempo in 1:steps
        config_new=aceptar(config_old,beta,n,m)
        out[tiempo+1]=energia_total(config_new,n,m)
        config_old,config_new=config_new,config_old
    end
    out
end

function magnetizaciones_t(beta,n::Int,m::Int,steps)
    out=zeros(steps+1)
    config_old=conf_aleatoria(n,m)
    out[1]=magnetizacion(config_old)
    for tiempo in 1:steps
        config_new=aceptar(config_old,beta,n,m)
        out[tiempo+1]=magnetizacion(config_new)
        config_old,config_new=config_new,config_old
    end
    out
end

function energiapromedio_beta(n::Int64,m::Int64)
    out=zeros(51)
    for i in 1:51
        out[i]=mean(energias_t((i-1)*0.1,n,m,100_000))
    end
    out
end

function magnetizacionpromedio_beta(n::Int64,m::Int64)
    out=zeros(51)
    for i in 1:51
        out[i]=mean(magnetizaciones_t((i-1)*0.1,n,m,100_000))
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
