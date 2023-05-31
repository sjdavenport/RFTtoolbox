function map = omega_index(n)
  map = [1,2,4];
  for i = 2:n
    map = [map, 4*(i-1)+[1,2,4]];
  end
end