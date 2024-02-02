   use LA_CONSTANTS, &
   only: wp=>sp, zero=>szero, half=>shalf, one=>sone, &
         safmin=>ssafmin, safmax=>ssafmax
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp) :: c, f, g, r, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, fs, g1, gs, u, rtmin, rtmax
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, sign, sqrt
!  ..
!  .. Constants ..
   rtmin = sqrt( safmin )
   rtmax = sqrt( safmax/2 )
!  ..
!  .. Executable Statements ..
!
   f1 = abs( f )
   g1 = abs( g )
   if( g == zero ) then
      c = one
      s = zero
      r = f
   else if( f == zero ) then
      c = zero
      s = sign( one, g )
      r = g1
   else if( f1 > rtmin .and. f1 < rtmax .and. &
            g1 > rtmin .and. g1 < rtmax ) then
      d = sqrt( f*f + g*g )
      c = f1 / d
      r = sign( d, f )
      s = g / r
   else
      u = min( safmax, max( safmin, f1, g1 ) )
      fs = f / u
      gs = g / u
      d = sqrt( fs*fs + gs*gs )
      c = abs( fs ) / d
      r = sign( d, f )
      s = gs / r
      r = r*u
   end if
   return
end subroutine