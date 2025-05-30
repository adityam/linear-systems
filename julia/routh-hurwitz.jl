using Polynomials 

using Hyperscript, Symbolics, Printf, LinearAlgebra, Latexify

function RouthArray(polynomial; details=true, normalize=true, s1zero=true)
  @variables ϵ 
  coeffs = length(polynomial) 
  n_rows = coeffs
  n_cols = [ ceil(Int, idx/2) for idx=n_rows:-1:1 ]

  math(x) = @sprintf("\$%s\$", x)

  RA = zeros(Num, n_rows, n_cols[1])
  factors = zeros(Int, n_rows)
  # flag = 1 means row of zero
  # flag = 2 means previous row zero
  flags = zeros(Int, n_rows)

  function normalize_row(row)
    row_values = Symbolics.value.(RA[row, :])
    if normalize && n_cols[row] != 1 && all(typeof.(row_values) .== Int)
       factor = gcd(row_values)
       if factor != 1
          factors[row] = factor
          RA[row, :] = div.(row_values, factor)
       end
    end
  end
  

  ## Fill the first two rows
  for idx in 1:coeffs
    RA[2 - rem(idx,2), ceil(Int, idx/2)] = polynomial[idx]
  end
  normalize_row.([1,2])

  for row in 3:n_rows
   ## Handle row of zero
   check_row = RA[row-1,:] .== 0
   check_first = RA[row-1,1] == 0
   flag = (row != n_rows) || s1zero
   if (all(typeof.(check_row) .== Bool) && all(check_row) && flag)
     flags[row-1] = 2
     row_order  = coeffs - (row - 2)
     row_powers = row_order:-2:0
     RA[row-1, :] = RA[row-2, :] .* row_powers
     normalize_row(row-1)
   elseif (typeof(check_first) == Bool) && check_first
     flags[row-1] = 1
     RA[row - 1, 1] = ϵ
   end

   ## Normal case
   for col = 1:n_cols[row]
      D = [ RA[row-2, 1] RA[row-2, col+1]; RA[row-1, 1] RA[row-1, col+1] ]
      RA[row, col] = -det(D)//RA[row - 1, 1] |> expand |> simplify
   end

  end
  
  ## Print the table
  table = m("table")
  tr = m("tr")
  td = m("td")

  function normalized_value(row, col)
    result = "" 
    if factors[row] != 0
       result= @sprintf("\\cancel{%s} %s",
                        ltx(RA[row, col]*factors[row]),
                        ltx(RA[row, col]))
        
    else
       result = ltx(RA[row, col])
    end
  end

  function value(row, col)
    result = ""
    if flags[row] == 1 && col == 1
       result = @sprintf("\\cancel{0} %s", normalized_value(row,col))
    elseif flags[row] == 2
       result = @sprintf("\\cancel{0} %s", normalized_value(row,col))
    else
       result = normalized_value(row,col)
    end

    if details && row > 2
      @sprintf("\$\\displaystyle -\\frac{\\DET{ %s & %s \\\\ %s & %s}}{%s} = %s\$", 
               ltx(RA[row - 2, 1]), ltx(RA[row - 2, col+1]),
               ltx(RA[row - 1, 1]), ltx(RA[row - 1, col+1]),
               ltx(RA[row - 1, 1]), result )
    else
      math(result)
    end
  end

  function ltx(x)
    latexify(x, env=:raw)
  end


  html = table.routharray(
    [ tr( td(@sprintf("\$s^{%s}\$", coeffs - row )),
          [td(value(row, col)) for col = 1:n_cols[row]]) for row = 1:n_rows ])
  
  return (html, RA)
end

function sign_changes(RA)
  sign(x) = (x > 0) ? "\$+\$" : "\$-\$"

  first_col = Symbolics.value.(RA[:, 1])
  n_rows = size(RA,1)

  ## Print the table
  table = m("table")
  tr = m("tr")
  td = m("td")
  thead = m("thead")
  @variables ϵ 

  if any(SymbolicUtils.symtype.(first_col) .== Real) # zero in first column
    table.routhsigns(
      thead( td("Term"), td("\$ ε \\to 0^{+}\$"), td("\$ ε \\to 0^{-}\$" )),
    [ tr( td(@sprintf("\$s^{%s}\$", n_rows - row )),
          td(sign(substitute(RA[row, 1], Dict([ ϵ =>  1e-10])))),
          td(sign(substitute(RA[row, 1], Dict([ ϵ => -1e-10]))))
        ) for row = 1:n_rows ])
  else
    table.routhsigns(
      thead( td("Term"), td("Sign") ),
    [ tr( td(@sprintf("\$s^{%s}\$", n_rows - row )),
          td(sign(RA[row, 1]))
        ) for row = 1:n_rows ])
  end
    
end

