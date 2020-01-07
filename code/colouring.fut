

let A:[][]f32 = [[1.0, 1.0, 0.0, 1.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                 [0.0, 1.0, 0.0, 0.0, 1.0, 1.0],
                 [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 1.0, 1.0, 0.0]]



let T:[][]f32 = [[1.000000f32, 1.000000f32, 0.000000f32, 0.000000f32, 0.000000f32],
                 [1.000000f32, 1.000000f32, 1.000000f32, 0.000000f32, 0.000000f32],
                 [0.000000f32, 1.000000f32, 1.000000f32, 1.000000f32, 0.000000f32],
                 [0.000000f32, 0.000000f32, 1.000000f32, 1.000000f32, 1.000000f32],
                 [0.000000f32, 0.000000f32, 0.000000f32, 1.000000f32, 1.000000f32]]

  -- tabulate_2d 5 5 (\i j -> if i==j then 1 else 0)

let matrix_to_bipartite [m][n] (A:[m][n]f32) : [m][n]i32 =
  tabulate_2d m n (\i j -> if f32.abs A[i,j] < 0.000001 then 0 else 1)

let argmin [n] (v_i:i32) (X:[n]i32) : i32 =
  if n < 1
  then
    -1
  else
    let XI = zip (iota n) X
    let XI = filter (\(_,x) -> x != v_i) XI
    in reduce (\(i1, x1) (i2, x2) -> if x1 <= x2 then (i1, x1) else (i2, x2)) XI[0] XI[1:length XI] |> (.1)

let max_colors:i32 = 10
-- G represents an adjency matrix of a bipartite graph G = (V_1, V_2, E)
-- where vertexes of V_1 is rows and
-- vertex of V_2 are the columns
-- There is an edge if element (i,j) is 1
let greedy_distance_2_coloring [m][n] (G:[m][n]i32) =
  let color = replicate n 0
  let tmp = tabulate (n*m) (\_ -> 0)
  let forbiddenColors = replicate max_colors (n+1)

  let (color, _, _) = loop (color, forbiddenColors, tmp) for j < n do
  let N_v = filter (>=0) <| tabulate m (\i -> if G[i,j] == 1 then i else (-1))  --- Neighbors of v_i
  let colored_index_N_w = flatten <| map (\j' -> tabulate m (\i -> if tmp[j' * m + i] == 1 then (i, j') else (-1, j'))) N_v
  let colored_index_N_w = filter (\(i,_) -> i >= 0) colored_index_N_w -- Neighbors of w that are colored
  let color_idx = map (.1) colored_index_N_w |> map (\i -> color[i])  -- Color indexes of colored Neighbors of w
  let forbiddenColors' = scatter forbiddenColors color_idx (replicate (length color_idx) j)
  let idx = map (\i -> i * m + j) N_v
  let color[j] = argmin j forbiddenColors'
  let new_tmp = scatter tmp idx (replicate (length idx) 1i32)
  in (color, forbiddenColors', new_tmp)
  in color


-- 1. Convert into bipartite graph
-- 2. Apply coloring on Jacobian
let main =
  let x = matrix_to_bipartite T
  in greedy_distance_2_coloring x
