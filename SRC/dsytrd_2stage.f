      SUBROUTINE DSYTRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU, HOUS2, LHOUS2, WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT, UPLO;
      int                N, LDA, LWORK, LHOUS2, INFO;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      double             A( LDA, * ), TAU( * ), HOUS2( * ), WORK( * );
      // ..

*  =====================================================================
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTQ;
      int                KD, IB, LWMIN, LHMIN, LWRK, LDAB, WPOS, ABPOS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DSYTRD_SY2SB, DSYTRD_SB2ST
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      // EXTERNAL LSAME, ILAENV2STAGE
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO   = 0
      WANTQ  = LSAME( VECT, 'V' )
      UPPER  = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 ) .OR. ( LHOUS2 == -1 )

      // Determine the block size, the workspace size and the hous size.

      KD     = ILAENV2STAGE( 1, 'DSYTRD_2STAGE', VECT, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'DSYTRD_2STAGE', VECT, N, KD, -1, -1 )
      if ( N == 0 ) {
         LHMIN = 1
         LWMIN = 1
      } else {
         LHMIN = ILAENV2STAGE( 3, 'DSYTRD_2STAGE', VECT, N, KD, IB, -1 )
         LWMIN = ILAENV2STAGE( 4, 'DSYTRD_2STAGE', VECT, N, KD, IB, -1 )
      }

      if ( .NOT.LSAME( VECT, 'N' ) ) {
         INFO = -1
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LHOUS2.LT.LHMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO == 0 ) {
         HOUS2( 1 ) = LHMIN
         WORK( 1 )  = LWMIN
      }

      if ( INFO != 0 ) {
         xerbla('DSYTRD_2STAGE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      // Determine pointer position

      LDAB  = KD+1
      LWRK  = LWORK-LDAB*N
      ABPOS = 1
      WPOS  = ABPOS + LDAB*N
      dsytrd_sy2sb(UPLO, N, KD, A, LDA, WORK( ABPOS ), LDAB, TAU, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('DSYTRD_SY2SB', -INFO );
         RETURN
      }
      dsytrd_sb2st('Y', VECT, UPLO, N, KD, WORK( ABPOS ), LDAB, D, E, HOUS2, LHOUS2, WORK( WPOS ), LWRK, INFO );
      if ( INFO != 0 ) {
         xerbla('DSYTRD_SB2ST', -INFO );
         RETURN
      }


      WORK( 1 ) = LWMIN
      RETURN

      // End of DSYTRD_2STAGE

      }
