      void zlarot(LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT, XRIGHT ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LLEFT, LRIGHT, LROWS;
      int                LDA, NL;
      Complex         C, S, XLEFT, XRIGHT;
      // ..
      // .. Array Arguments ..
      Complex         A( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                IINC, INEXT, IX, IY, IYT, J, NT;
      Complex         TEMPX;
      // ..
      // .. Local Arrays ..
      Complex         XT( 2 ), YT( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      // .. Executable Statements ..

      // Set up indices, arrays for ends

      if ( LROWS ) {
         IINC = LDA;
         INEXT = 1;
      } else {
         IINC = 1;
         INEXT = LDA;
      }

      if ( LLEFT ) {
         NT = 1;
         IX = 1 + IINC;
         IY = 2 + LDA;
         XT[1] = A( 1 );
         YT[1] = XLEFT;
      } else {
         NT = 0;
         IX = 1;
         IY = 1 + INEXT;
      }

      if ( LRIGHT ) {
         IYT = 1 + INEXT + ( NL-1 )*IINC;
         NT = NT + 1;
         XT[NT] = XRIGHT;
         YT[NT] = A( IYT );
      }

      // Check for errors

      if ( NL < NT ) {
         xerbla('ZLAROT', 4 );
         return;
      }
      if ( LDA <= 0 || ( !LROWS && LDA < NL-NT ) ) {
         xerbla('ZLAROT', 8 );
         return;
      }

      // Rotate

      // ZROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S

      for (J = 0; J <= NL - NT - 1; J++) { // 10
         TEMPX = C*A( IX+J*IINC ) + S*A( IY+J*IINC );
         A[IY+J*IINC] = -DCONJG( S )*A( IX+J*IINC ) + DCONJG( C )*A( IY+J*IINC );
         A[IX+J*IINC] = TEMPX;
      } // 10

      // ZROT( NT, XT,1, YT,1, C, S ) with complex C, S

      for (J = 1; J <= NT; J++) { // 20
         TEMPX = C*XT( J ) + S*YT( J );
         YT[J] = -DCONJG( S )*XT( J ) + DCONJG( C )*YT( J );
         XT[J] = TEMPX;
      } // 20

      // Stuff values back into XLEFT, XRIGHT, etc.

      if ( LLEFT ) {
         A[1] = XT( 1 );
         XLEFT = YT( 1 );
      }

      if ( LRIGHT ) {
         XRIGHT = XT( NT );
         A[IYT] = YT( NT );
      }

      return;
      }