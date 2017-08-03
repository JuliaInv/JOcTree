
function getF( x, y, z )
p2 = 2*pi
f = sin.(p2*x) .* sin.(p2*y) .* sin.(p2*z)
return f
end  # function getF


# Derivatives

function getdFdX( x, y, z )
p2 = 2*pi
dfx = (p2*cos.(p2*x)) .* sin.(p2*y) .* sin.(p2*z)
return dfx
end  # function getdFdX


function getdFdY( x, y, z )
p2 = 2*pi
dfy = sin.(p2*x) .* (p2*cos.(p2*y)) .* sin.(p2*z)
return dfy
end  # function getdFdY


function getdFdZ( x, y, z )
p2 = 2*pi
dfz = sin.(p2*x) .* sin.(p2*y) .* (p2*cos.(p2*z))
return dfz
end  # function getdFdZ
