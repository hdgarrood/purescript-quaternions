module Test.Main where

import Prelude
import Data.Monoid
import Data.Array as Array
import Data.Foldable (foldMap, and)
import Data.Vector as Vec
import Data.Vector3 (Vec3)
import Data.Vector3 as Vec3
import Test.QuickCheck
import Test.QuickCheck.Gen as Gen
import Test.QuickCheck.Arbitrary
import Control.Monad.Eff.Console (log)

import Data.Quaternion hiding (approxEq)
import Data.Quaternion as Quaternion
import Data.Quaternion.Rotation (Rotation)
import Data.Quaternion.Rotation as Rotation

newtype ArbQ = ArbQ (Quaternion Number)

-- A generator of Number values between -100 and 100.
smallNum = map (\x -> (x * 200.0) - 100.0) arbitrary

instance arbQ :: Arbitrary ArbQ where
  arbitrary =
    (\a b c d -> ArbQ (Quaternion a b c d))
      <$> smallNum
      <*> smallNum
      <*> smallNum
      <*> smallNum

newtype ArbRot = ArbRot (Rotation Number)

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

epsilon = 0.00000001

approxEq x y = Math.abs (x - y) < epsilon

vApproxEq x y =
  case Vec.toArray x, Vec.toArray y of
    xs, ys ->
      and (Array.zipWith approxEq xs ys)

qApproxEq = Quaternion.approxEq epsilon

newtype LargeArray a = LargeArray (Array a)

instance arbLargeArray :: Arbitrary a => Arbitrary (LargeArray a) where
  arbitrary =
    map LargeArray (Gen.chooseInt 100 1000 >>= \x -> Gen.vectorOf x arbitrary)

main = do
  let showR = Rotation.showAngleAxis

  log "fromAngleAxis and toAngleAxis are approximate inverses"
  quickCheck \(ArbRot p) ->
    let
      q = Rotation.fromAngleAxis (Rotation.toAngleAxis p)
    in
      qApproxEq (Rotation.toQuaternion p) (Rotation.toQuaternion q)
      <?> ("p: " <> show p <> ", q: " <> show q)

  log "Rotations are versors"
  quickCheck \(ArbRot p) ->
    let
      n = norm (Rotation.toQuaternion p)
    in
      approxEq n 1.0 <?> ("p: " <> show p <> ", n: " <> show n)

  log "Numerical stability"
  quickCheck \(LargeArray xs) ->
    let
      product = foldMap runArbRot xs
      n = norm (Rotation.toQuaternion product)
    in
      approxEq n 1.0 <?> ("count: " <> show (Array.length xs) <> ", n: " <> show n)

  log "act is a group action: identity"
  quickCheck \(ArbV3 v) ->
    Rotation.act mempty v == v <?> ("v: " <> show v)

  log "act is a group action: compatibility"
  quickCheck \(ArbV3 v) (ArbRot p) (ArbRot q) ->
    vApproxEq (Rotation.act (p <> q) v) (Rotation.act p (Rotation.act q v))
    <?> ("p: " <> showR p <> ", q: " <> showR q <> ", v: " <> show v)
