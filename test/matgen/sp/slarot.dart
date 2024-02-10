      void slarot(LROWS, LLEFT, LRIGHT, NL, C, S, final Matrix<double> A, final int LDA, XLEFT, final int XRIGHT) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               LLEFT, LRIGHT, LROWS;
      int                LDA, NL;
      double               C, S, XLEFT, XRIGHT;
      double               A( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                IINC, INEXT, IX, IY, IYT, NT;
      double               XT( 2 ), YT( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SROT, XERBLA

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
         xerbla('SLAROT', 4 );
         return;
      }
      if ( LDA <= 0 || ( !LROWS && LDA < NL-NT ) ) {
         xerbla('SLAROT', 8 );
         return;
      }

      // Rotate

      srot(NL-NT, A( IX ), IINC, A( IY ), IINC, C, S );
      srot(NT, XT, 1, YT, 1, C, S );

      // Stuff values back into XLEFT, XRIGHT, etc.

      if ( LLEFT ) {
         A[1] = XT( 1 );
         XLEFT = YT( 1 );
      }

      if ( LRIGHT ) {
         XRIGHT = XT( NT );
         A[IYT] = YT( NT );
      }

      }
