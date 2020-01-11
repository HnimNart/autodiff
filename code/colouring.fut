import "ad"

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

let min (x1:i32) (x2:i32) = if x1 < x2 then x1 else x2

let max_colors:i32 = 10
-- G represents an adjency matrix of a bipartite graph G = (V_1, V_2, E)
-- where vertexes of V_1 is rows and
-- vertex of V_2 are the columns
-- There is an edge if element (i,j) is 1
-- Implementation is sequential in each coloring
let greedy_distance_2_coloring [m][n] (G:[m][n]i32) =
  let coloring = replicate n 0
  let colored = replicate (n*m) (false) -- to keep track of Neighbors that has been coloured
  let forbiddenColors = replicate (min max_colors n) (n+1)
  let (coloring, _, _) =
    loop (coloring, forbiddenColors, colored) for j < n do
    let N_v_index = filter (>=0) <| tabulate m (\i -> if G[i,j] == 1 then i else (-1))  --- Neighbors of v_i, i.e. w's
    let colored_index_N_w = flatten <| map (\j' -> tabulate m (\i -> if unsafe colored[j' * m + i] then i else (-1))) N_v_index
    let colored_index_N_w = filter (>= 0) colored_index_N_w -- Neighbors of w that are colored
    let color_idx = map (\i -> coloring[i]) colored_index_N_w   -- Color indexes of colored Neighbors of w
    let forbiddenColors' = scatter forbiddenColors color_idx (replicate (length color_idx) j)
    let idx = map (\i -> i * m + j) N_v_index
    let coloring[j] = argmin j forbiddenColors'
    let colored' = scatter colored idx (replicate (length idx) true)
    in (coloring, forbiddenColors', colored')
  in coloring



module f32_dual = mk_dual f32


let f xs =
  let n = length xs
  in map (\i -> if i == 0 then f32_dual.(xs[0] + xs[1])
                else if i == n - 1 then
                     let x1 = xs[n - 2]
                     let x2 = xs[n - 1]
                     in f32_dual.(x1 + x2)
                else
                     let x1 = xs[i-1]
                     let x2 = xs[i]
                     let x3 = xs[i+1]
                     in f32_dual.(x1 + x2 + x3)) (iota (n))


let maximum [n] (xs:[n]i32): i32 =
  reduce (\x1 x2 -> if x1 < x2 then x2 else x1) (-1) xs

let compute_seed_matrix [n] (coloring:[n]i32) :[][]f32 =
  let max_col = maximum coloring
  in tabulate max_col (\i -> map (\c -> if c == i then 1 else 0) coloring)


-- 1. Convert into bipartite graph
-- 2. Apply coloring on Jacobian
let main =
  let dim = 5
  let input = tabulate dim (\i -> f32.i32 i)
  let dual_input = tabulate dim (\i -> map2 f32_dual.make_dual input (tabulate dim (\j -> f32.bool(j==i))))
  let res = map f dual_input
  let jacobian = map (map (.2)) res
  let G = matrix_to_bipartite jacobian
  let coloring = greedy_distance_2_coloring G
  let compressed_seed = compute_seed_matrix coloring
  let dual_input = map (\row -> map2 f32_dual.make_dual input row) compressed_seed
  let res = map f dual_input

  let compressed_J =  res |> map unzip |> unzip |> (.2)
  let J = tabulate dim (\i -> -- restore uncompressed Jacobian
                          let col = compressed_J[coloring[i]] in
                          let g_col = G[i] in map2 (\c mask -> c * f32.i32 mask) col g_col)
  in J
