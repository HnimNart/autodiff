import "ad"
import "lib/github.com/athas/vector/vspace"


let diag (X:[]f32) =
  let len = length X
  let elem = len ** 2
  let index  = map (\i -> i * len + i) (0..<len)
  let retval = scatter (replicate elem (0)) index X
  in unflatten len len retval


module fs (M: real) = {
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

  let project cam_r cam_c cam_f cam_p cam_k X =
    let Xcam = rodrigues_rotate_point
               (point_3d cam_r)
               (X v3d.- (point_3d cam_c))
    let distorted = radial_distort (point_2d cam_k)
                                   (p2e Xcam)
    in point_2d cam_p v2d.+ v2d.scale cam_f distorted

  let compute_zach_weight_error w =
    M.(i32 1 - w*w)

  let compute_reproj_err cam_r cam_c cam_f cam_p cam_k X_points w feat =
    v2d.scale w <| project cam_r cam_c cam_f cam_p cam_k X_points v2d.- feat

}

module f32_dual = mk_dual f32
module fs_f32_dual = fs f32_dual

let cam_n = 11
let param_len = cam_n + 4


-- cams: 11 cameras in format [r1 r2 r3 C1 C2 C3 f u0 v0 k1 k2] OK
--       r1, r2, r3 are angle - axis rotation parameters(Rodrigues)
--       [C1 C2 C3]' is the camera center
--       f is the focal length in pixels
--       [u0 v0]' is the principal point
--       k1, k2 are radial distortion parameters
-- X: 3*m points
-- obs: 2*p observations (pairs cameraIdx, pointIdx)
-- feats: 2*p features (x,y coordinates corresponding to observations)
-- reproj_err: 2*p errors of observations
let compute_reproj_err_wrapper (target:fs_f32_dual.point_2d)
                               (input:[param_len]f32_dual.t) =
  let cam_r = input[:3]
  let cam_c = input[3:6]
  let cam_f = input[6]
  let cam_p = input[7:9]
  let cam_k = input[9:11]
  let X = {x=input[11], y=input[12], z=input[13]}
  let w = input[14]
  in fs_f32_dual.compute_reproj_err cam_r cam_c cam_f cam_p cam_k X w target

-- Computes reprojection error for a single data-point X and camera
let compute_reproj_err (cam: [cam_n]f32) (X: [3]f32) (w:f32) (feat:[2]f32) =
  let input = cam ++ X ++ [w] -- Strect input into one long vector
  let seed_matrix = diag (replicate param_len 1)
  let input_dual = map (map2 f32_dual.make_dual input) seed_matrix
  -- Target does not need derivatives
  let target_dual = fs_f32_dual.point_2d (map f32_dual.inject feat)

  in map (compute_reproj_err_wrapper target_dual) input_dual


let main [n][m][p] (cams: [n][cam_n]f32)
                   (X: [m][3]f32)
                   (w: [p]f32)
                   (obs: [p][2]i32)
                   (feats: [p][2]f32) =
  let cam_idx = obs[:, 0]
  let point_idx = obs[:, 1]
  let res = compute_reproj_err cams[0] X[0] w[0] feats[0]
  let w_err = map
       -- (map (\i -> unsafe cams[i]) obs[:,0])
       -- (map (\i -> unsafe X[i]) obs[:,1])
       -- w feats

-- x1 * x2 + sin(x1)
-- let main =
--   let x1 = f32_dual.make_dual 3 1
--   let x2 = f32_dual.make_dual 2 0
--   let sinx1 = f32_dual.sin x1
--   in ((f32_dual.(x1 * x2 + sinx1)).2, (f32_dual.(x1 * x2)).2)
