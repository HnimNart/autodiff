import "ad"
import "lib/github.com/athas/vector/vspace"


module M = f32
module f32_dual = mk_dual M

module v3d = mk_vspace_3d M
type point_3d = v3d.vector
let point_3d a: point_3d = {x=a[0], y=a[1], z=a[2]}

module v2d = mk_vspace_2d M
type point_2d = v2d.vector
let point_2d a: point_2d = {x=a[0], y=a[1]}

let rodrigues_rotate_point (rot: point_3d) (X: point_3d) =
  let sqtheta = v3d.quadrance rot in
  if M.(sqtheta != i32 0) then
  let theta = M.sqrt sqtheta
  let costheta = M.cos theta
  let sintheta = M.sin theta
  let theta_inv = M.(i32 1 / theta)

  let w = v3d.scale theta_inv rot
  let w_cross_X = v3d.cross w X
  let tmp = M.(v3d.dot w X * (i32 1 - costheta))

  in v3d.(scale costheta X +
          scale sintheta w_cross_X +
          scale tmp w)
  else v3d.(X + cross rot X)


let radial_distort (rad_params: point_2d) (proj: point_2d) =
  let rsq = v2d.quadrance proj
  let L = M.(i32 1 + rad_params.x * rsq + rad_params.y * rsq * rsq)
  in v2d.scale L proj


let p2e (X:point_3d) : point_2d =
  v2d.scale M.(i32 1/X.z) {x=X.x, y=X.y}




let main =
  let x1 = f32_dual.make_dual 2 1
  let x2 = f32_dual.make_dual 3 0
  let x1x2 = f32_dual.(x1 * x2)
  let sinx1 = f32_dual.sin x1
  in f32_dual.(x1x2 + sinx1)
