      SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANS;
      int                IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, LWORK, M, N;
      double             DIF, SCALE;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * ), WORK( * )
      // ..
*
*  =====================================================================
*  Replaced various illegal calls to CCOPY by calls to CLASET.
*  Sven Hammarling, 1/5/02.
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = (0.0D+0, 0.0D+0) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOTRAN;
      int                I, IE, IFUNC, IROUND, IS, ISOLVE, J, JE, JS, K, LINFO, LWMIN, MB, NB, P, PQ, Q;
      double             DSCALE, DSUM, SCALE2, SCALOC;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZLACPY, ZLASET, ZSCAL, ZTGSY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, SQRT
      // ..
      // .. Executable Statements ..
*
      // Decode and test input parameters
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( NOTRAN ) THEN
         IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.4 ) ) THEN
            INFO = -2
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( M.LE.0 ) THEN
            INFO = -3
         ELSE IF( N.LE.0 ) THEN
            INFO = -4
         ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -6
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -8
         ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
         ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
            INFO = -12
         ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
            INFO = -16
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( NOTRAN ) THEN
            IF( IJOB.EQ.1 .OR. IJOB.EQ.2 ) THEN
               LWMIN = MAX( 1, 2*M*N )
            ELSE
               LWMIN = 1
            END IF
         ELSE
            LWMIN = 1
         END IF
         WORK( 1 ) = LWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -20
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSYL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         SCALE = 1
         IF( NOTRAN ) THEN
            IF( IJOB.NE.0 ) THEN
               DIF = 0
            END IF
         END IF
         RETURN
      END IF
*
      // Determine  optimal block sizes MB and NB
*
      MB = ILAENV( 2, 'ZTGSYL', TRANS, M, N, -1, -1 )
      NB = ILAENV( 5, 'ZTGSYL', TRANS, M, N, -1, -1 )
*
      ISOLVE = 1
      IFUNC = 0
      IF( NOTRAN ) THEN
         IF( IJOB.GE.3 ) THEN
            IFUNC = IJOB - 2
            CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
            CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
         ELSE IF( IJOB.GE.1 .AND. NOTRAN ) THEN
            ISOLVE = 2
         END IF
      END IF
*
      IF( ( MB.LE.1 .AND. NB.LE.1 ) .OR. ( MB.GE.M .AND. NB.GE.N ) ) THEN
*
         // Use unblocked Level 2 solver
*
         DO 30 IROUND = 1, ISOLVE
*
            SCALE = ONE
            DSCALE = ZERO
            DSUM = ONE
            PQ = M*N
            CALL ZTGSY2( TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DSUM, DSCALE, INFO )
            IF( DSCALE.NE.ZERO ) THEN
               IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                  DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               ELSE
                  DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               END IF
            END IF
            IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
               IF( NOTRAN ) THEN
                  IFUNC = IJOB
               END IF
               SCALE2 = SCALE
               CALL ZLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL ZLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
            ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
               CALL ZLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL ZLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            END IF
   30    CONTINUE
*
         RETURN
*
      END IF
*
      // Determine block structure of A
*
      P = 0
      I = 1
   40 CONTINUE
      IF( I.GT.M ) GO TO 50
      P = P + 1
      IWORK( P ) = I
      I = I + MB
      IF( I.GE.M ) GO TO 50
      GO TO 40
   50 CONTINUE
      IWORK( P+1 ) = M + 1
      IF( IWORK( P ).EQ.IWORK( P+1 ) ) P = P - 1
*
      // Determine block structure of B
*
      Q = P + 1
      J = 1
   60 CONTINUE
      IF( J.GT.N ) GO TO 70
*
      Q = Q + 1
      IWORK( Q ) = J
      J = J + NB
      IF( J.GE.N ) GO TO 70
      GO TO 60
*
   70 CONTINUE
      IWORK( Q+1 ) = N + 1
      IF( IWORK( Q ).EQ.IWORK( Q+1 ) ) Q = Q - 1
*
      IF( NOTRAN ) THEN
         DO 150 IROUND = 1, ISOLVE
*
            // Solve (I, J) - subsystem
                // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
                // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
            // for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
*
            PQ = 0
            SCALE = ONE
            DSCALE = ZERO
            DSUM = ONE
            DO 130 J = P + 2, Q
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               DO 120 I = P, 1, -1
                  IS = IWORK( I )
                  IE = IWORK( I+1 ) - 1
                  MB = IE - IS + 1
                  CALL ZTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, LINFO )
                  IF( LINFO.GT.0 ) INFO = LINFO
                  PQ = PQ + MB*NB
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 K = 1, JS - 1
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                         CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
   80                CONTINUE
                     DO 90 K = JS, JE
                        CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                         CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
   90                CONTINUE
                     DO 100 K = JS, JE
                        CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ), C( IE+1, K ), 1 )                         CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ), F( IE+1, K ), 1 )
  100                CONTINUE
                     DO 110 K = JE + 1, N
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                         CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
  110                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
*
                  // Substitute R(I,J) and L(I,J) into remaining equation.
*
                  IF( I.GT.1 ) THEN
                     CALL ZGEMM( 'N', 'N', IS-1, NB, MB, DCMPLX( -ONE, ZERO ), A( 1, IS ), LDA, C( IS, JS ), LDC, DCMPLX( ONE, ZERO ), C( 1, JS ), LDC )                      CALL ZGEMM( 'N', 'N', IS-1, NB, MB, DCMPLX( -ONE, ZERO ), D( 1, IS ), LDD, C( IS, JS ), LDC, DCMPLX( ONE, ZERO ), F( 1, JS ), LDF )
                  END IF
                  IF( J.LT.Q ) THEN
                     CALL ZGEMM( 'N', 'N', MB, N-JE, NB, DCMPLX( ONE, ZERO ), F( IS, JS ), LDF, B( JS, JE+1 ), LDB, DCMPLX( ONE, ZERO ), C( IS, JE+1 ), LDC )                      CALL ZGEMM( 'N', 'N', MB, N-JE, NB, DCMPLX( ONE, ZERO ), F( IS, JS ), LDF, E( JS, JE+1 ), LDE, DCMPLX( ONE, ZERO ), F( IS, JE+1 ), LDF )
                  END IF
  120          CONTINUE
  130       CONTINUE
            IF( DSCALE.NE.ZERO ) THEN
               IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                  DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               ELSE
                  DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               END IF
            END IF
            IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
               IF( NOTRAN ) THEN
                  IFUNC = IJOB
               END IF
               SCALE2 = SCALE
               CALL ZLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL ZLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
            ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
               CALL ZLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL ZLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            END IF
  150    CONTINUE
      ELSE
*
         // Solve transposed (I, J)-subsystem
             // A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J)
             // R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J)
         // for I = 1,2,..., P; J = Q, Q-1,..., 1
*
         SCALE = ONE
         DO 210 I = 1, P
            IS = IWORK( I )
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 200 J = Q, P + 2, -1
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               CALL ZTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, LINFO )
               IF( LINFO.GT.0 ) INFO = LINFO
               IF( SCALOC.NE.ONE ) THEN
                  DO 160 K = 1, JS - 1
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                      CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
  160             CONTINUE
                  DO 170 K = JS, JE
                     CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                      CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
  170             CONTINUE
                  DO 180 K = JS, JE
                     CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ), C( IE+1, K ), 1 )                      CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ), F( IE+1, K ), 1 )
  180             CONTINUE
                  DO 190 K = JE + 1, N
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ), 1 )                      CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ), 1 )
  190             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
*
               // Substitute R(I,J) and L(I,J) into remaining equation.
*
               IF( J.GT.P+2 ) THEN
                  CALL ZGEMM( 'N', 'C', MB, JS-1, NB, DCMPLX( ONE, ZERO ), C( IS, JS ), LDC, B( 1, JS ), LDB, DCMPLX( ONE, ZERO ), F( IS, 1 ), LDF )                   CALL ZGEMM( 'N', 'C', MB, JS-1, NB, DCMPLX( ONE, ZERO ), F( IS, JS ), LDF, E( 1, JS ), LDE, DCMPLX( ONE, ZERO ), F( IS, 1 ), LDF )
               END IF
               IF( I.LT.P ) THEN
                  CALL ZGEMM( 'C', 'N', M-IE, NB, MB, DCMPLX( -ONE, ZERO ), A( IS, IE+1 ), LDA, C( IS, JS ), LDC, DCMPLX( ONE, ZERO ), C( IE+1, JS ), LDC )                   CALL ZGEMM( 'C', 'N', M-IE, NB, MB, DCMPLX( -ONE, ZERO ), D( IS, IE+1 ), LDD, F( IS, JS ), LDF, DCMPLX( ONE, ZERO ), C( IE+1, JS ), LDC )
               END IF
  200       CONTINUE
  210    CONTINUE
      END IF
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
      // End of ZTGSYL
*
      END
