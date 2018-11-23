module Test.Main where

import Prelude

import Data.Array as Array
import Data.Foldable (foldl, foldr, foldMap, sum)
import Data.Quaternion (Quaternion(..))
import Data.Quaternion as Quaternion
import Data.Quaternion.Rotation (Rotation)
import Data.Quaternion.Rotation as Rotation
import Data.Quaternion.Vec3 (Vec3)
import Data.Quaternion.Vec3 as Vec3
import Effect (Effect)
import Effect.Console (log)
import Math as Math
import Partial.Unsafe (unsafePartial)
import Test.QuickCheck (quickCheck, (<?>))
import Test.QuickCheck.Arbitrary (class Arbitrary, arbitrary)
import Test.QuickCheck.Gen as Gen

-- Have a 3x3 matrix (in column-major order) act on a vector in 3D space via
-- multiplication
mat3mul :: Partial => Array Number -> Vec3 Number -> Vec3 Number
mat3mul [a11, a21, a31, a12, a22, a32, a13, a23, a33] v =
  Vec3.vec3
    (Vec3.dot (Vec3.vec3 a11 a12 a13) v)
    (Vec3.dot (Vec3.vec3 a21 a22 a23) v)
    (Vec3.dot (Vec3.vec3 a31 a32 a33) v)

newtype ArbQ = ArbQ (Quaternion Number)

-- A generator of Number values between -100 and 100.
smallNum :: Gen.Gen Number
smallNum = map (\x -> (x * 200.0) - 100.0) arbitrary

instance arbQ :: Arbitrary ArbQ where
  arbitrary =
    (\a b c d -> ArbQ (Quaternion a b c d))
      <$> smallNum
      <*> smallNum
      <*> smallNum
      <*> smallNum

newtype ArbRot = ArbRot Rotation

runArbRot :: ArbRot -> Rotation
runArbRot (ArbRot p) = p

instance arbRot :: Arbitrary ArbRot where
  arbitrary = go <$> smallNum <*> arbitrary
    where
    go angle (ArbV3 axis) =
      ArbRot (Rotation.fromAngleAxis { angle, axis })

newtype ArbV3 = ArbV3 (Vec3 Number)

instance arbV3 :: Arbitrary ArbV3 where
  arbitrary =
    (\x y z -> ArbV3 (Vec3.vec3 x y z))
      <$> smallNum
      <*> smallNum
      <*> smallNum

epsilon :: Number
epsilon = 0.00000001

approxEq :: Number -> Number -> Boolean
approxEq x y = Math.abs (x - y) < epsilon

vApproxEq :: Vec3 Number -> Vec3 Number -> Boolean
vApproxEq x y =
  Vec3.magnitude (Vec3.vsub x y) < epsilon

qApproxEq :: Quaternion Number -> Quaternion Number -> Boolean
qApproxEq = Quaternion.approxEq epsilon

-- Approximate equality for rotations
rApproxEq :: Rotation -> Rotation -> Boolean
rApproxEq p q =
  let
    p' = Rotation.toQuaternion p
    q' = Rotation.toQuaternion q
  in
    -- note that negating a rotation quaternion actually gives you another
    -- quaternion representing the exact same rotation
    (qApproxEq p' q' || qApproxEq (negate p') q')

newtype LargeArray a = LargeArray (Array a)

instance arbLargeArray :: Arbitrary a => Arbitrary (LargeArray a) where
  arbitrary =
    map LargeArray (Gen.chooseInt 100 1000 >>= \x -> Gen.vectorOf x arbitrary)

matDistance :: Array Number -> Array Number -> Number
matDistance x y =
  let
    magnitude x' = sum (Array.zipWith (*) x' x')
  in
    magnitude (Array.zipWith (-) x y)

main :: Effect Unit
main = do
  log "Foldable instance agrees with arrays"
  quickCheck \(ArbQ p@(Quaternion a b c d)) ->
    -- Use a non-associative operator here so that foldl and foldr are
    -- distinguishable
    let
      p' = [a, b, c, d]
    in
      (foldl (-) 0.0 p == foldl (-) 0.0 p')
      && (foldr (-) 0.0 p == foldr (-) 0.0 p')
      <?> ("p: " <> show p)

  let showR = Rotation.showAngleAxis
  log "fromAngleAxis and toAngleAxis are approximate inverses"
  quickCheck \(ArbRot p) ->
    let
      q = Rotation.fromAngleAxis (Rotation.toAngleAxis p)
    in
      rApproxEq p q
      <?> ("p: " <> show p <> ", q: " <> show q)

  log "Rotation matrices agree with rotations"
  quickCheck \(ArbRot p) (ArbV3 v) ->
    let
      w1 = Rotation.act p v
      w2 = unsafePartial (mat3mul (Rotation.toRotationMatrix p) v)
    in
      vApproxEq w1 w2
      <?> ("p: " <> show p <> ", v: " <> show v)

  log "fromRotationMatrix and toRotationMatrix are approximate inverses"
  quickCheck \(ArbRot p) ->
    let
      q = unsafePartial (Rotation.fromRotationMatrix (Rotation.toRotationMatrix p))
    in
      rApproxEq p q
      <?> ("p: " <> show p <> ", q: " <> show q)

  log "Rotations are versors"
  quickCheck \(ArbRot p) ->
    let
      n = Quaternion.norm (Rotation.toQuaternion p)
    in
      approxEq n 1.0 <?> ("p: " <> show p <> ", n: " <> show n)

  log "Numerical stability"
  quickCheck \(LargeArray xs) ->
    let
      product = foldMap runArbRot xs
      n = Quaternion.norm (Rotation.toQuaternion product)
    in
      approxEq n 1.0 <?> ("count: " <> show (Array.length xs) <> ", n: " <> show n)

  log "act is a group action: identity"
  quickCheck \(ArbV3 v) ->
    Rotation.act mempty v == v <?> ("v: " <> show v)

  log "act is a group action: compatibility"
  quickCheck \(ArbV3 v) (ArbRot p) (ArbRot q) ->
    vApproxEq (Rotation.act (p <> q) v) (Rotation.act p (Rotation.act q v))
    <?> ("p: " <> showR p <> ", q: " <> showR q <> ", v: " <> show v)

  -- log "(toMat4 <<< fromAngleAxis) equivalent to Data.Matrix.makeRotate"
  -- quickCheck \(ArbV3 axis) angle ->
  --   let
  --     dist = matDistance
  --             (Mat4.makeRotate angle axis)
  --             (Rotation.toMat4 (Rotation.fromAngleAxis { angle, axis }))
  --   in
  --     approxEq dist 0.0 <?> ("angle: " <> show angle <> ", axis: " <> show axis)
