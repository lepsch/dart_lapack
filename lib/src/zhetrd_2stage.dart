      void zhetrd_2stage(final int VECT, final int UPLO, final int N, final Matrix<double> A, final int LDA, final int D, final int E, final int TAU, final int HOUS2, final int LHOUS2, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             VECT, UPLO;
      int                N, LDA, LWORK, LHOUS2, INFO;
      double             D( * ), E( * );
      Complex         A( LDA, * ), TAU( * ), HOUS2( * ), WORK( * );
      // ..

// =====================================================================
      bool               LQUERY, UPPER, WANTQ;
      int                KD, IB, LWMIN, LHMIN, LWRK, LDAB, WPOS, ABPOS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHETRD_HE2HB, ZHETRD_HB2ST
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      // EXTERNAL lsame, ILAENV2STAGE

      // Test the input parameters

      INFO   = 0;
      WANTQ  = lsame( VECT, 'V' );
      UPPER  = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 ) || ( LHOUS2 == -1 );

      // Determine the block size, the workspace size and the hous size.

      KD     = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', VECT, N, -1, -1, -1 );
      IB     = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', VECT, N, KD, -1, -1 );
      if ( N == 0 ) {
         LHMIN = 1;
         LWMIN = 1;
      } else {
         LHMIN = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', VECT, N, KD, IB, -1 );
         LWMIN = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', VECT, N, KD, IB, -1 );
      }

      if ( !lsame( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LHOUS2 < LHMIN && !LQUERY ) {
         INFO = -10;
      } else if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO == 0 ) {
         HOUS2[1] = LHMIN;
         WORK[1] = LWMIN;
      }

      if ( INFO != 0 ) {
         xerbla('ZHETRD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK[1] = 1;
         return;
      }

      // Determine pointer position

      LDAB  = KD+1;
      LWRK  = LWORK-LDAB*N;
      ABPOS = 1;
      WPOS  = ABPOS + LDAB*N;
      zhetrd_he2hb(UPLO, N, KD, A, LDA, WORK( ABPOS ), LDAB, TAU, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('ZHETRD_HE2HB', -INFO );
         return;
      }
      zhetrd_hb2st('Y', VECT, UPLO, N, KD, WORK( ABPOS ), LDAB, D, E, HOUS2, LHOUS2, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('ZHETRD_HB2ST', -INFO );
         return;
      }


      WORK[1] = LWMIN;
      }
