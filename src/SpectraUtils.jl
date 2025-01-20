module SpectraUtils
using Unitful
import PhysicalConstants.CODATA2018: c_0, h, k_B

export wavenumber, wavelength, gigahertz, terahertz, electronvolt
export dopplerwidth

Wavenumber = Union{Quantity{T, Unitful.𝐋^-1, U}, Level{L, S, Quantity{T, Unitful.𝐋^-1, U}} where {L, S}} where {T, U}

wavenumber(λ::Unitful.Length) = 1/λ |> u"cm^-1"
wavenumber(λ::Unitful.Length, Δλ::Unitful.Length) = (1(λ), 1/(λ) - 1/(λ+Δλ)) |> u"cm^-1"
wavenumber(E::Unitful.Energy) = E / (h*c_0) |> u"cm^-1"
wavenumber(E::Unitful.Energy, ΔE::Unitful.Energy) = (E/(h*c_0), ΔE/(h*c_0)) |> u"cm^-1"
wavenumber(f::Unitful.Frequency) = f / (c_0) |> u"cm^-1"
wavenumber(f::Unitful.Frequency, Δf::Unitful.Frequency) = (f/(c_0), Δf/(c_0)) |> u"cm^-1"

wavelength(ν::Wavenumber) = 1/ν |> u"nm"
wavelength(ν::Wavenumber, Δν::Wavenumber) = 1/ν |> u"nm", 1/ν - 1/(ν+Δν) |> u"nm"
wavelength(x...) = wavelength(wavenumber(x...))

gigahertz(ν::Wavenumber) = c_0 * ν |> u"GHz"
gigahertz(ν::Wavenumber, Δν::Wavenumber) = c_0*ν |> u"GHz", c_0*Δν |> u"GHz"
gigahertz(x...) = gigahertz(wavenumber(x...))

terahertz(x...) = gigahertz(x...) |> u"THz"

electronvolt(ν::Wavenumber) = h*c_0 * ν |> u"eV"
electronvolt(ν::Wavenumber, Δν::Wavenumber) = h*c_0*ν |> u"eV", h*c_0*Δν |> u"eV"
electronvolt(x...) = electronvolt(wavenumber(x...))


dopplerwidth(ν::Wavenumber, m::Unitful.Mass, T::Unitful.Temperature) = sqrt(8*k_B * T * log(2)/ m/c_0^2) * ν
dopplerwidth(x::Unitful.Energy, args...) = electronvolt(dopplerwidth(wavenumber(x), args...) )
dopplerwidth(x::Unitful.Frequency, args...) = gigahertz(dopplerwidth(wavenumber(x), args...) )
dopplerwidth(x::Unitful.Length, args...) = wavelength(wavenumber(x), dopplerwidth(wavenumber(x), args...))[2]
end
