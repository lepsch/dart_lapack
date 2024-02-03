      void chetrd_2stage(VECT, UPLO, N, A, LDA, D, E, TAU, HOUS2, LHOUS2, WORK, LWORK, INFO ) {

      // IMPLICIT NONE

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT, UPLO;
      int                N, LDA, LWORK, LHOUS2, INFO;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * );
      COMPLEX            A( LDA, * ), TAU( * ), HOUS2( * ), WORK( * );
      // ..

// =====================================================================
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTQ;
      int                KD, IB, LWMIN, LHMIN, LWRK, LDAB, WPOS, ABPOS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHETRD_HE2HB, CHETRD_HB2ST
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SROUNDUP_LWORK;
      // EXTERNAL LSAME, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO   = 0;
      WANTQ  = LSAME( VECT, 'V' );
      UPPER  = LSAME( UPLO, 'U' );
      LQUERY = ( LWORK == -1 ) || ( LHOUS2 == -1 );

      // Determine the block size, the workspace size and the hous size.

      KD     = ILAENV2STAGE( 1, 'CHETRD_2STAGE', VECT, N, -1, -1, -1 );
      IB     = ILAENV2STAGE( 2, 'CHETRD_2STAGE', VECT, N, KD, -1, -1 );
      if ( N == 0 ) {
         LHMIN = 1;
         LWMIN = 1;
      } else {
         LHMIN = ILAENV2STAGE( 3, 'CHETRD_2STAGE', VECT, N, KD, IB, -1 );
         LWMIN = ILAENV2STAGE( 4, 'CHETRD_2STAGE', VECT, N, KD, IB, -1 );
      }

      if ( !LSAME( VECT, 'N' ) ) {
         INFO = -1;
      } else if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
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
         HOUS2( 1 ) = SROUNDUP_LWORK( LHMIN );
         WORK( 1 )  = SROUNDUP_LWORK( LWMIN );
      }

      if ( INFO != 0 ) {
         xerbla('CHETRD_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK( 1 ) = 1;
         return;
      }

      // Determine pointer position

      LDAB  = KD+1;
      LWRK  = LWORK-LDAB*N;
      ABPOS = 1;
      WPOS  = ABPOS + LDAB*N;
      chetrd_he2hb(UPLO, N, KD, A, LDA, WORK( ABPOS ), LDAB, TAU, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('CHETRD_HE2HB', -INFO );
         return;
      }
      chetrd_hb2st('Y', VECT, UPLO, N, KD, WORK( ABPOS ), LDAB, D, E, HOUS2, LHOUS2, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('CHETRD_HB2ST', -INFO );
         return;
      }


      WORK( 1 ) = SROUNDUP_LWORK( LWMIN );
      return;

      // End of CHETRD_2STAGE

      }
