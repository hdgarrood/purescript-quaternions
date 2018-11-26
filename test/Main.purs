module Test.Main where

import Prelude

import Data.Array as Array
import Data.Foldable (class Foldable, foldl, foldr, foldMap, sum)
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
import Test.QuickCheck as QC
import Test.QuickCheck.Arbitrary (class Arbitrary, arbitrary)
import Test.QuickCheck.Gen as Gen

foldResults :: forall f. Foldable f => f QC.Result -> QC.Result
foldResults = foldr combine QC.Success
  where
  combine QC.Success r = r
  combine r@(QC.Failed str) _ = r

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

newtype SmallNum = SmallNum Number

instance arbitrarySmallNum :: Arbitrary SmallNum where
  arbitrary = SmallNum <$> smallNum

-- A generator of Number values between 0 and 1/100000.
verySmallNum :: Gen.Gen Number
verySmallNum = map (_ / 100000.0) arbitrary

newtype VerySmallNum = VerySmallNum Number

instance arbitraryVerySmallNum :: Arbitrary VerySmallNum where
  arbitrary = VerySmallNum <$> verySmallNum

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
      ArbRot (Rotation.fromAxisAngle { angle, axis })

newtype ArbV3 = ArbV3 (Vec3 Number)

instance arbV3 :: Arbitrary ArbV3 where
  arbitrary =
    (\x y z -> ArbV3 (Vec3.vec3 x y z))
      <$> smallNum
      <*> smallNum
      <*> smallNum

epsilon :: Number
epsilon = 1e-8

approxEq :: Number -> Number -> Boolean
approxEq x y = Math.abs (x - y) < epsilon

vApproxEq :: Vec3 Number -> Vec3 Number -> Boolean
vApproxEq x y =
  Vec3.norm (Vec3.vsub x y) < epsilon

qApproxEq :: Quaternion Number -> Quaternion Number -> Boolean
qApproxEq p q =
  let
    -- Scale the tolerance based on the norms of the arguments
    tolerance = Math.max (Quaternion.norm p) (Quaternion.norm q) / 1e12
  in
    Quaternion.approxEq tolerance p q

-- Approximate equality for rotations
rApproxEq :: Rotation -> Rotation -> Boolean
rApproxEq = Rotation.approxEq epsilon

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
  log "fromReal/fromVector identity"
  quickCheck \(ArbQ q) ->
    let
      q' = Quaternion.fromReal (Quaternion.realPart q) +
            Quaternion.fromVector (Quaternion.vectorPart q)
    in
      q == q'
      <?> show { q, q' }

  log "scalarMul agrees with Quaternion multiplication"
  quickCheck \(ArbQ p) k ->
    Quaternion.scalarMul k p == (Quaternion k 0.0 0.0 0.0) * p
    <?> show { p, k }

  log "versor/scalarMul identity"
  quickCheck \(ArbQ q) ->
    let
      q' = Quaternion.scalarMul (Quaternion.norm q) (Quaternion.versor q)
    in
      qApproxEq q q'
      <?> show { q, q' }

  log "Foldable instance agrees with arrays"
  quickCheck \(ArbQ p@(Quaternion a b c d)) ->
    -- Use a non-associative operator here so that foldl and foldr are
    -- distinguishable
    let
      p' = [a, b, c, d]
    in
      (foldl (-) 0.0 p == foldl (-) 0.0 p')
      && (foldr (-) 0.0 p == foldr (-) 0.0 p')
      <?> show { p }

  log "Quaternion exponential agrees with real exponential"
  quickCheck \x ->
    Quaternion.exp (Quaternion.fromReal x) == Quaternion (Math.exp x) 0.0 0.0 0.0
    <?> show { x }

  log "exp(p+q) = exp(p)*exp(q) when p, q commute"
  quickCheck \x1 y1 x2 y2 ->
    let
      p = Quaternion x1 0.0 y1 0.0
      q = Quaternion x2 0.0 y2 0.0
    in
      qApproxEq (Quaternion.exp (p + q)) (Quaternion.exp p * Quaternion.exp q)
      <?> show { p, q }

  log "exp <<< log == identity (approximately)"
  quickCheck \(ArbQ q) ->
    let
      expLogQ = Quaternion.exp (Quaternion.log q)
    in
      qApproxEq expLogQ q
      <?> show { q, expLogQ }

  log "exp agrees with the complex exponential function"
  quickCheck \(SmallNum x') (SmallNum y) ->
    let
      x = x' / 5.0 -- this gives us an x in [-20, 20]
      w = Quaternion x y 0.0 0.0
      z@(Quaternion _ _ e1 e2) = Quaternion.exp w
      z' = Quaternion.scalarMul (Math.exp x)
            (Quaternion.fromReal (Math.cos y) +
              Quaternion.i * Quaternion.fromReal (Math.sin y))
    in
      foldResults
        [ e1 == 0.0 && e2 == 0.0 <?> show { w, z, msg: "expected z to be complex" }
        , qApproxEq z z' <?> show { w, z, z' }
        ]

  log "log agrees with the principal value of the complex logarithm"
  quickCheck \(SmallNum x) (SmallNum y) ->
    let
      w = Quaternion x y 0.0 0.0
      z@(Quaternion w' x' y' z') = Quaternion.log w
    in
      foldResults
        [ y' == 0.0 && z' == 0.0 <?> show { w, z, msg: "expected z to be complex" }
        , -Math.pi <= x' && x' <= Math.pi <?> show { w, z, msg: "expected x' to be in [-pi, pi]" }
        ]

  log "fromAxisAngle and toAxisAngle are approximate inverses"
  quickCheck \(ArbRot p) ->
    let
      q = Rotation.fromAxisAngle (Rotation.toAxisAngle p)
    in
      rApproxEq p q
      <?> show { p, q }

  log "toAxisAngle returns a unit-length axis"
  quickCheck \(ArbRot p) ->
    let
      { axis } = Rotation.toAxisAngle p
    in
      approxEq (Vec3.norm axis) 1.0
      <?> show { p }

  log "fromAxisAngle doesn't mind about the axis' magnitude"
  quickCheck \(ArbV3 v) theta k ->
    let
      p = Rotation.fromAxisAngle { angle: theta, axis: v }
      q = Rotation.fromAxisAngle { angle: theta, axis: Vec3.scalarMul k v }
    in
      rApproxEq p q
      <?> show { v, k }

  log "Rotation matrices agree with rotations"
  quickCheck \(ArbRot p) (ArbV3 v) ->
    let
      w1 = Rotation.act p v
      w2 = unsafePartial (mat3mul (Rotation.toRotationMatrix p) v)
    in
      vApproxEq w1 w2
      <?> show { p, v }

  log "fromRotationMatrix and toRotationMatrix are approximate inverses"
  quickCheck \(ArbRot p) ->
    let
      q = unsafePartial (Rotation.fromRotationMatrix (Rotation.toRotationMatrix p))
    in
      rApproxEq p q
      <?> show { p, q }

  log "fromRotationMatrix avoids instability issues when the quaternion has one very small component"
  quickCheck \(VerySmallNum e) a b c ->
    let
      testWith :: Quaternion Number -> QC.Result
      testWith p' =
        let
          p = Rotation.fromQuaternion p'
          q = unsafePartial (Rotation.fromRotationMatrix (Rotation.toRotationMatrix p))
        in
          rApproxEq p q
          <?> show { p, q }

      z = 0.0
    in
      foldResults
        [ testWith (Quaternion e a b c)
        , testWith (Quaternion z a b c)
        , testWith (Quaternion a e b c)
        , testWith (Quaternion a z b c)
        , testWith (Quaternion a b e c)
        , testWith (Quaternion a b z c)
        , testWith (Quaternion a b e c)
        , testWith (Quaternion a b e z)
        ]

  log "Rotations are versors"
  quickCheck \(ArbRot p) ->
    let
      n = Quaternion.norm (Rotation.toQuaternion p)
    in
      approxEq n 1.0 <?> show { p, n }

  log "Numerical stability"
  quickCheck \(LargeArray xs) ->
    let
      product = foldMap runArbRot xs
      n = Quaternion.norm (Rotation.toQuaternion product)
    in
      approxEq n 1.0 <?> show { count: Array.length xs, n }

  log "act is a group action: identity"
  quickCheck \(ArbV3 v) ->
    Rotation.act mempty v == v <?> show { v }

  log "act is a group action: compatibility"
  quickCheck \(ArbV3 v) (ArbRot p) (ArbRot q) ->
    vApproxEq (Rotation.act (p <> q) v) (Rotation.act p (Rotation.act q v))
    <?> show { p, q, v }
  -- log "(toMat4 <<< fromAxisAngle) equivalent to Data.Matrix.makeRotate"
  -- quickCheck \(ArbV3 axis) angle ->
  --   let
  --     dist = matDistance
  --             (Mat4.makeRotate angle axis)
  --             (Rotation.toMat4 (Rotation.fromAxisAngle { angle, axis }))
  --   in
  --     approxEq dist 0.0 <?> ("angle: " <> show angle <> ", axis: " <> show axis)
