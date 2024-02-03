      SUBROUTINE ZHETRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU, HOUS2, LHOUS2, WORK, LWORK, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             VECT, UPLO;
      int                N, LDA, LWORK, LHOUS2, INFO;
*     ..
*     .. Array Arguments ..
      double             D( * ), E( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), HOUS2( * ), WORK( * )
*     ..
*
*  =====================================================================
*     ..
*     .. Local Scalars ..
      bool               LQUERY, UPPER, WANTQ;
      int                KD, IB, LWMIN, LHMIN, LWRK, LDAB, WPOS, ABPOS;
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA, ZHETRD_HE2HB, ZHETRD_HB2ST
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      // EXTERNAL LSAME, ILAENV2STAGE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO   = 0
      WANTQ  = LSAME( VECT, 'V' )
      UPPER  = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 ) .OR. ( LHOUS2.EQ.-1 )
*
*     Determine the block size, the workspace size and the hous size.
*
      KD     = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', VECT, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', VECT, N, KD, -1, -1 )
      IF( N.EQ.0 ) THEN
         LHMIN = 1
         LWMIN = 1
      ELSE
         LHMIN = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', VECT, N, KD, IB, -1 )
         LWMIN = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', VECT, N, KD, IB, -1 )
      END IF
*
      IF( .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LHOUS2.LT.LHMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         HOUS2( 1 ) = LHMIN
         WORK( 1 )  = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD_2STAGE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine pointer position
*
      LDAB  = KD+1
      LWRK  = LWORK-LDAB*N
      ABPOS = 1
      WPOS  = ABPOS + LDAB*N
      CALL ZHETRD_HE2HB( UPLO, N, KD, A, LDA, WORK( ABPOS ), LDAB, TAU, WORK( WPOS ), LWRK, INFO )
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD_HE2HB', -INFO )
         RETURN
      END IF
      CALL ZHETRD_HB2ST( 'Y', VECT, UPLO, N, KD, WORK( ABPOS ), LDAB, D, E, HOUS2, LHOUS2, WORK( WPOS ), LWRK, INFO )
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD_HB2ST', -INFO )
         RETURN
      END IF
*
*
      WORK( 1 )  = LWMIN
      RETURN
*
*     End of ZHETRD_2STAGE
*
      END
