      SUBROUTINE DLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT, XRIGHT )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LLEFT, LRIGHT, LROWS;
      int                LDA, NL;
      double             C, S, XLEFT, XRIGHT;
      // ..
      // .. Array Arguments ..
      double             A( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                IINC, INEXT, IX, IY, IYT, NT;
      // ..
      // .. Local Arrays ..
      double             XT( 2 ), YT( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DROT, XERBLA
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

      if ( NL.LT.NT ) {
         CALL XERBLA( 'DLAROT', 4 )
         RETURN
      }
      if ( LDA.LE.0 .OR. ( .NOT.LROWS .AND. LDA.LT.NL-NT ) ) {
         CALL XERBLA( 'DLAROT', 8 )
         RETURN
      }

      // Rotate

      CALL DROT( NL-NT, A( IX ), IINC, A( IY ), IINC, C, S )
      CALL DROT( NT, XT, 1, YT, 1, C, S )

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

      // End of DLAROT

      }
