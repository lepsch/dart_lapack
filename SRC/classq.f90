   use LA_CONSTANTS, &
      only: wp=>sp, zero=>szero, one=>sone, &
            sbig=>ssbig, ssml=>sssml, tbig=>stbig, tsml=>stsml
   use LA_XISNAN
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scale, sumsq
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, ymax, ymin
!  ..
!
!  Quick return if possible
!
   if( LA_ISNAN(scale) .or. LA_ISNAN(sumsq) ) return
   if( sumsq == zero ) scale = one
   if( scale == zero ) then
      scale = one
      sumsq = zero
   end if
   if (n <= 0) then
      return
   end if
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(real(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ax = abs(aimag(x(ix)))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Put the existing sum of squares into one of the accumulators
!
   if( sumsq > zero ) then
      ax = scale*sqrt( sumsq )
      if (ax > tbig) then
         if (scale > one) then
            scale = scale * sbig
            abig = abig + scale * (scale * sumsq)
         else
            ! sumsq > tbig^2 => (sbig * (sbig * sumsq)) is representable
            abig = abig + scale * (scale * (sbig * (sbig * sumsq)))
         end if
      else if (ax < tsml) then
         if (notbig) then
            if (scale < one) then
               scale = scale * ssml
               asml = asml + scale * (scale * sumsq)
            else
               ! sumsq < tsml^2 => (ssml * (ssml * sumsq)) is representable
               asml = asml + scale * (scale * (ssml * (ssml * sumsq)))
            end if
         end if
      else
         amed = amed + scale * (scale * sumsq)
      end if
   end if
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if (amed > zero .or. LA_ISNAN(amed)) then
         abig = abig + (amed*sbig)*sbig
      end if
      scale = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if (amed > zero .or. LA_ISNAN(amed)) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scale = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scale = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range or zero
!
      scale = one
      sumsq = amed
   end if
   return
end subroutine