      SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL, LRWMIN, LWMIN, M, SMLSIZ, START
      REAL               EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANST, SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANST, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CLACPY, CLACRM, CLAED0, CSTEQR, CSWAP, SLASCL, SLASET, SSTEDC, SSTEQR, SSTERF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         SMLSIZ = ILAENV( 9, 'CSTEDC', ' ', 0, 0, 0, 0 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 2*( N - 1 )
         ELSE IF( ICOMPZ.EQ.1 ) THEN
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWMIN = N*N
            LRWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1
            LRWMIN = 1 + 4*N + 2*N**2
            LIWMIN = 3 + 5*N
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures SSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default. If the conditional clause is removed, then
*     information on the size of workspace needs to be changed.
*
*     If COMPZ = 'N', use SSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SSTERF( N, D, E, INFO )
         GO TO 70
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
*
         CALL CSTEQR( COMPZ, N, D, E, Z, LDZ, RWORK, INFO )
*
      ELSE
*
*        If COMPZ = 'I', we simply call SSTEDC instead.
*
         IF( ICOMPZ.EQ.2 ) THEN
            CALL SLASET( 'Full', N, N, ZERO, ONE, RWORK, N )
            LL = N*N + 1
            CALL SSTEDC( 'I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1, IWORK, LIWORK, INFO )
            DO 20 J = 1, N
               DO 10 I = 1, N
                  Z( I, J ) = RWORK( ( J-1 )*N+I )
   10          CONTINUE
   20       CONTINUE
            GO TO 70
         END IF
*
*        From now on, only option left to be handled is COMPZ = 'V',
*        i.e. ICOMPZ = 1.
*
*        Scale.
*
         ORGNRM = SLANST( 'M', N, D, E )
         IF( ORGNRM.EQ.ZERO ) GO TO 70
*
         EPS = SLAMCH( 'Epsilon' )
*
         START = 1
*
*        while ( START <= N )
*
   30    CONTINUE
         IF( START.LE.N ) THEN
*
*           Let FINISH be the position of the next subdiagonal entry
*           such that E( FINISH ) <= TINY or FINISH = N if no such
*           subdiagonal exists.  The matrix identified by the elements
*           between START and FINISH constitutes an independent
*           sub-problem.
*
            FINISH = START
   40       CONTINUE
            IF( FINISH.LT.N ) THEN
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 40
               END IF
            END IF
*
*           (Sub) Problem determined.  Compute its size and solve it.
*
            M = FINISH - START + 1
            IF( M.GT.SMLSIZ ) THEN
*
*              Scale.
*
               ORGNRM = SLANST( 'M', M, D( START ), E( START ) )
               CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO )                CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO )
*
               CALL CLAED0( N, M, D( START ), E( START ), Z( 1, START ), LDZ, WORK, N, RWORK, IWORK, INFO )
               IF( INFO.GT.0 ) THEN
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 70
               END IF
*
*              Scale back.
*
               CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO )
*
            ELSE
               CALL SSTEQR( 'I', M, D( START ), E( START ), RWORK, M, RWORK( M*M+1 ), INFO )                CALL CLACRM( N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N, RWORK( M*M+1 ) )
               CALL CLACPY( 'A', N, M, WORK, N, Z( 1, START ), LDZ )
               IF( INFO.GT.0 ) THEN
                  INFO = START*( N+1 ) + FINISH
                  GO TO 70
               END IF
            END IF
*
            START = FINISH + 1
            GO TO 30
         END IF
*
*        endwhile
*
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 60 II = 2, N
           I = II - 1
           K = I
           P = D( I )
           DO 50 J = II, N
              IF( D( J ).LT.P ) THEN
                 K = J
                 P = D( J )
              END IF
   50      CONTINUE
           IF( K.NE.I ) THEN
              D( K ) = D( I )
              D( I ) = P
              CALL CSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
           END IF
   60    CONTINUE
      END IF
*
   70 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CSTEDC
*
      END
