      SUBROUTINE CHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      int                K, LDA, N
      String             TRANS, TRANSR, UPLO;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      REAL               ONE, ZERO
      COMPLEX            CZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, NORMALTRANSR, NISODD, NOTRANS
      int                INFO, NROWA, J, NK, N1, N2
      COMPLEX            CALPHA, CBETA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CHERK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, CMPLX
*     ..
*     .. Executable Statements ..
*
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )
*
      IF( NOTRANS ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
*
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHFRK ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
*     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
*     done (it is in CHERK for example) and left in the general case.
*
      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND. ( BETA.EQ.ONE ) ) )RETURN
*
      IF( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ZERO ) ) THEN
         DO J = 1, ( ( N*( N+1 ) ) / 2 )
            C( J ) = CZERO
         END DO
         RETURN
      END IF
*
      CALPHA = CMPLX( ALPHA, ZERO )
      CBETA = CMPLX( BETA, ZERO )
*
*     C is N-by-N.
*     If N is odd, set NISODD = .TRUE., and N1 and N2.
*     If N is even, NISODD = .FALSE., and NK.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         NISODD = .FALSE.
         NK = N / 2
      ELSE
         NISODD = .TRUE.
         IF( LOWER ) THEN
            N2 = N / 2
            N1 = N - N2
         ELSE
            N1 = N / 2
            N2 = N - N1
         END IF
      END IF
*
      IF( NISODD ) THEN
*
*        N is odd
*
         IF( NORMALTRANSR ) THEN
*
*           N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              N is odd, TRANSR = 'N', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
*
                  CALL CHERK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL CHERK( 'U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N )                   CALL CGEMM( 'N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N )
*
               ELSE
*
*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
*
                  CALL CHERK( 'L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL CHERK( 'U', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N )                   CALL CGEMM( 'C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( N1+1 ), N )
*
               END IF
*
            ELSE
*
*              N is odd, TRANSR = 'N', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
*
                  CALL CHERK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL CHERK( 'U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N )                   CALL CGEMM( 'N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, CBETA, C( 1 ), N )
*
               ELSE
*
*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
*
                  CALL CHERK( 'L', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL CHERK( 'U', 'C', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N )                   CALL CGEMM( 'C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, CBETA, C( 1 ), N )
*
               END IF
*
            END IF
*
         ELSE
*
*           N is odd, and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              N is odd, TRANSR = 'C', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
*
                  CALL CHERK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL CHERK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 )                   CALL CGEMM( 'N', 'C', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, CBETA, C( N1*N1+1 ), N1 )
*
               ELSE
*
*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
*
                  CALL CHERK( 'U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL CHERK( 'L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 )                   CALL CGEMM( 'C', 'N', N1, N2, K, CALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, CBETA, C( N1*N1+1 ), N1 )
*
               END IF
*
            ELSE
*
*              N is odd, TRANSR = 'C', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
*
                  CALL CHERK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL CHERK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL CGEMM( 'N', 'C', N2, N1, K, CALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 )
*
               ELSE
*
*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
*
                  CALL CHERK( 'U', 'C', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL CHERK( 'L', 'C', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL CGEMM( 'C', 'N', N2, N1, K, CALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), N2 )
*
               END IF
*
            END IF
*
         END IF
*
      ELSE
*
*        N is even
*
         IF( NORMALTRANSR ) THEN
*
*           N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              N is even, TRANSR = 'N', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
*
                  CALL CHERK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL CHERK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 )                   CALL CGEMM( 'N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 )
*
               ELSE
*
*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
*
                  CALL CHERK( 'L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL CHERK( 'U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 )                   CALL CGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( NK+2 ), N+1 )
*
               END IF
*
            ELSE
*
*              N is even, TRANSR = 'N', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
*
                  CALL CHERK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL CHERK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL CGEMM( 'N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( 1 ), N+1 )
*
               ELSE
*
*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
*
                  CALL CHERK( 'L', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL CHERK( 'U', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL CGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( 1 ), N+1 )
*
               END IF
*
            END IF
*
         ELSE
*
*           N is even, and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              N is even, TRANSR = 'C', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
*
                  CALL CHERK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL CHERK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK )                   CALL CGEMM( 'N', 'C', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK )
*
               ELSE
*
*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
*
                  CALL CHERK( 'U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL CHERK( 'L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK )                   CALL CGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, CBETA, C( ( ( NK+1 )*NK )+1 ), NK )
*
               END IF
*
            ELSE
*
*              N is even, TRANSR = 'C', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
*
                  CALL CHERK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL CHERK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL CGEMM( 'N', 'C', NK, NK, K, CALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK )
*
               ELSE
*
*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
*
                  CALL CHERK( 'U', 'C', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL CHERK( 'L', 'C', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL CGEMM( 'C', 'N', NK, NK, K, CALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, CBETA, C( 1 ), NK )
*
               END IF
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of CHFRK
*
      END
