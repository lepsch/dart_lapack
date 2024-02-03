      SUBROUTINE CLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT, XRIGHT )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LLEFT, LRIGHT, LROWS;
      int                LDA, NL;
      COMPLEX            C, S, XLEFT, XRIGHT
      // ..
      // .. Array Arguments ..
      COMPLEX            A( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                IINC, INEXT, IX, IY, IYT, J, NT;
      COMPLEX            TEMPX
      // ..
      // .. Local Arrays ..
      COMPLEX            XT( 2 ), YT( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Set up indices, arrays for ends

      if ( LROWS ) {
         IINC = LDA
         INEXT = 1
      } else {
         IINC = 1
         INEXT = LDA
      }

      if ( LLEFT ) {
         NT = 1
         IX = 1 + IINC
         IY = 2 + LDA
         XT( 1 ) = A( 1 )
         YT( 1 ) = XLEFT
      } else {
         NT = 0
         IX = 1
         IY = 1 + INEXT
      }

      if ( LRIGHT ) {
         IYT = 1 + INEXT + ( NL-1 )*IINC
         NT = NT + 1
         XT( NT ) = XRIGHT
         YT( NT ) = A( IYT )
      }

      // Check for errors

      if ( NL < NT ) {
         xerbla('CLAROT', 4 );
         RETURN
      }
      if ( LDA <= 0 || ( !LROWS && LDA < NL-NT ) ) {
         xerbla('CLAROT', 8 );
         RETURN
      }

      // Rotate

      // CROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S

      for (J = 0; J <= NL - NT - 1; J++) { // 10
         TEMPX = C*A( IX+J*IINC ) + S*A( IY+J*IINC )
         A( IY+J*IINC ) = -CONJG( S )*A( IX+J*IINC ) + CONJG( C )*A( IY+J*IINC )
         A( IX+J*IINC ) = TEMPX
      } // 10

      // CROT( NT, XT,1, YT,1, C, S ) with complex C, S

      for (J = 1; J <= NT; J++) { // 20
         TEMPX = C*XT( J ) + S*YT( J )
         YT( J ) = -CONJG( S )*XT( J ) + CONJG( C )*YT( J )
         XT( J ) = TEMPX
      } // 20

      // Stuff values back into XLEFT, XRIGHT, etc.

      if ( LLEFT ) {
         A( 1 ) = XT( 1 )
         XLEFT = YT( 1 )
      }

      if ( LRIGHT ) {
         XRIGHT = XT( NT )
         A( IYT ) = YT( NT )
      }

      RETURN

      // End of CLAROT

      }
