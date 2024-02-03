      SUBROUTINE STGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, LWORK, M, N;
      REAL               DIF, SCALE
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * ), WORK( * )
      // ..

*  =====================================================================
*  Replaced various illegal calls to SCOPY by calls to SLASET.
*  Sven Hammarling, 1/5/02.

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOTRAN;
      int                I, IE, IFUNC, IROUND, IS, ISOLVE, J, JE, JS, K, LINFO, LWMIN, MB, NB, P, PPQQ, PQ, Q;
      REAL               DSCALE, DSUM, SCALE2, SCALOC
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASET, SSCAL, STGSY2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test input parameters

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -1
      } else if ( NOTRAN ) {
         if ( ( IJOB.LT.0 ) .OR. ( IJOB.GT.4 ) ) {
            INFO = -2
         }
      }
      if ( INFO.EQ.0 ) {
         if ( M.LE.0 ) {
            INFO = -3
         } else if ( N.LE.0 ) {
            INFO = -4
         } else if ( LDA.LT.MAX( 1, M ) ) {
            INFO = -6
         } else if ( LDB.LT.MAX( 1, N ) ) {
            INFO = -8
         } else if ( LDC.LT.MAX( 1, M ) ) {
            INFO = -10
         } else if ( LDD.LT.MAX( 1, M ) ) {
            INFO = -12
         } else if ( LDE.LT.MAX( 1, N ) ) {
            INFO = -14
         } else if ( LDF.LT.MAX( 1, M ) ) {
            INFO = -16
         }
      }

      if ( INFO.EQ.0 ) {
         if ( NOTRAN ) {
            if ( IJOB.EQ.1 .OR. IJOB.EQ.2 ) {
               LWMIN = MAX( 1, 2*M*N )
            } else {
               LWMIN = 1
            }
         } else {
            LWMIN = 1
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -20
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'STGSYL', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         SCALE = 1
         if ( NOTRAN ) {
            if ( IJOB.NE.0 ) {
               DIF = 0
            }
         }
         RETURN
      }

      // Determine optimal block sizes MB and NB

      MB = ILAENV( 2, 'STGSYL', TRANS, M, N, -1, -1 )
      NB = ILAENV( 5, 'STGSYL', TRANS, M, N, -1, -1 )

      ISOLVE = 1
      IFUNC = 0
      if ( NOTRAN ) {
         if ( IJOB.GE.3 ) {
            IFUNC = IJOB - 2
            CALL SLASET( 'F', M, N, ZERO, ZERO, C, LDC )
            CALL SLASET( 'F', M, N, ZERO, ZERO, F, LDF )
         } else if ( IJOB.GE.1 .AND. NOTRAN ) {
            ISOLVE = 2
         }
      }

      if ( ( MB.LE.1 .AND. NB.LE.1 ) .OR. ( MB.GE.M .AND. NB.GE.N ) ) {

         DO 30 IROUND = 1, ISOLVE

            // Use unblocked Level 2 solver

            DSCALE = ZERO
            DSUM = ONE
            PQ = 0
            CALL STGSY2( TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DSUM, DSCALE, IWORK, PQ, INFO )
            if ( DSCALE.NE.ZERO ) {
               if ( IJOB.EQ.1 .OR. IJOB.EQ.3 ) {
                  DIF = SQRT( REAL( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               } else {
                  DIF = SQRT( REAL( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               }
            }

            if ( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) {
               if ( NOTRAN ) {
                  IFUNC = IJOB
               }
               SCALE2 = SCALE
               CALL SLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL SLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL SLASET( 'F', M, N, ZERO, ZERO, C, LDC )
               CALL SLASET( 'F', M, N, ZERO, ZERO, F, LDF )
            } else if ( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) {
               CALL SLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL SLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            }
   30    CONTINUE

         RETURN
      }

      // Determine block structure of A

      P = 0
      I = 1
   40 CONTINUE
      IF( I.GT.M ) GO TO 50
      P = P + 1
      IWORK( P ) = I
      I = I + MB
      IF( I.GE.M ) GO TO 50       IF( A( I, I-1 ).NE.ZERO ) I = I + 1
      GO TO 40
   50 CONTINUE

      IWORK( P+1 ) = M + 1
      IF( IWORK( P ).EQ.IWORK( P+1 ) ) P = P - 1

      // Determine block structure of B

      Q = P + 1
      J = 1
   60 CONTINUE
      IF( J.GT.N ) GO TO 70
      Q = Q + 1
      IWORK( Q ) = J
      J = J + NB
      IF( J.GE.N ) GO TO 70       IF( B( J, J-1 ).NE.ZERO ) J = J + 1
      GO TO 60
   70 CONTINUE

      IWORK( Q+1 ) = N + 1
      IF( IWORK( Q ).EQ.IWORK( Q+1 ) ) Q = Q - 1

      if ( NOTRAN ) {

         DO 150 IROUND = 1, ISOLVE

            // Solve (I, J)-subsystem
                // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
                // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
            // for I = P, P - 1,..., 1; J = 1, 2,..., Q

            DSCALE = ZERO
            DSUM = ONE
            PQ = 0
            SCALE = ONE
            DO 130 J = P + 2, Q
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               DO 120 I = P, 1, -1
                  IS = IWORK( I )
                  IE = IWORK( I+1 ) - 1
                  MB = IE - IS + 1
                  PPQQ = 0
                  CALL STGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, IWORK( Q+2 ), PPQQ, LINFO )
                  IF( LINFO.GT.0 ) INFO = LINFO

                  PQ = PQ + PPQQ
                  if ( SCALOC.NE.ONE ) {
                     DO 80 K = 1, JS - 1
                        CALL SSCAL( M, SCALOC, C( 1, K ), 1 )
                        CALL SSCAL( M, SCALOC, F( 1, K ), 1 )
   80                CONTINUE
                     DO 90 K = JS, JE
                        CALL SSCAL( IS-1, SCALOC, C( 1, K ), 1 )
                        CALL SSCAL( IS-1, SCALOC, F( 1, K ), 1 )
   90                CONTINUE
                     DO 100 K = JS, JE
                        CALL SSCAL( M-IE, SCALOC, C( IE+1, K ), 1 )
                        CALL SSCAL( M-IE, SCALOC, F( IE+1, K ), 1 )
  100                CONTINUE
                     DO 110 K = JE + 1, N
                        CALL SSCAL( M, SCALOC, C( 1, K ), 1 )
                        CALL SSCAL( M, SCALOC, F( 1, K ), 1 )
  110                CONTINUE
                     SCALE = SCALE*SCALOC
                  }

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I.GT.1 ) {
                     CALL SGEMM( 'N', 'N', IS-1, NB, MB, -ONE, A( 1, IS ), LDA, C( IS, JS ), LDC, ONE, C( 1, JS ), LDC )                      CALL SGEMM( 'N', 'N', IS-1, NB, MB, -ONE, D( 1, IS ), LDD, C( IS, JS ), LDC, ONE, F( 1, JS ), LDF )
                  }
                  if ( J.LT.Q ) {
                     CALL SGEMM( 'N', 'N', MB, N-JE, NB, ONE, F( IS, JS ), LDF, B( JS, JE+1 ), LDB, ONE, C( IS, JE+1 ), LDC )                      CALL SGEMM( 'N', 'N', MB, N-JE, NB, ONE, F( IS, JS ), LDF, E( JS, JE+1 ), LDE, ONE, F( IS, JE+1 ), LDF )
                  }
  120          CONTINUE
  130       CONTINUE
            if ( DSCALE.NE.ZERO ) {
               if ( IJOB.EQ.1 .OR. IJOB.EQ.3 ) {
                  DIF = SQRT( REAL( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               } else {
                  DIF = SQRT( REAL( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               }
            }
            if ( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) {
               if ( NOTRAN ) {
                  IFUNC = IJOB
               }
               SCALE2 = SCALE
               CALL SLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL SLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL SLASET( 'F', M, N, ZERO, ZERO, C, LDC )
               CALL SLASET( 'F', M, N, ZERO, ZERO, F, LDF )
            } else if ( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) {
               CALL SLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL SLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            }
  150    CONTINUE

      } else {

         // Solve transposed (I, J)-subsystem
              // A(I, I)**T * R(I, J)  + D(I, I)**T * L(I, J)  =  C(I, J)
              // R(I, J)  * B(J, J)**T + L(I, J)  * E(J, J)**T = -F(I, J)
         // for I = 1,2,..., P; J = Q, Q-1,..., 1

         SCALE = ONE
         DO 210 I = 1, P
            IS = IWORK( I )
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 200 J = Q, P + 2, -1
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               CALL STGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, IWORK( Q+2 ), PPQQ, LINFO )
               IF( LINFO.GT.0 ) INFO = LINFO
               if ( SCALOC.NE.ONE ) {
                  DO 160 K = 1, JS - 1
                     CALL SSCAL( M, SCALOC, C( 1, K ), 1 )
                     CALL SSCAL( M, SCALOC, F( 1, K ), 1 )
  160             CONTINUE
                  DO 170 K = JS, JE
                     CALL SSCAL( IS-1, SCALOC, C( 1, K ), 1 )
                     CALL SSCAL( IS-1, SCALOC, F( 1, K ), 1 )
  170             CONTINUE
                  DO 180 K = JS, JE
                     CALL SSCAL( M-IE, SCALOC, C( IE+1, K ), 1 )
                     CALL SSCAL( M-IE, SCALOC, F( IE+1, K ), 1 )
  180             CONTINUE
                  DO 190 K = JE + 1, N
                     CALL SSCAL( M, SCALOC, C( 1, K ), 1 )
                     CALL SSCAL( M, SCALOC, F( 1, K ), 1 )
  190             CONTINUE
                  SCALE = SCALE*SCALOC
               }

               // Substitute R(I, J) and L(I, J) into remaining equation.

               if ( J.GT.P+2 ) {
                  CALL SGEMM( 'N', 'T', MB, JS-1, NB, ONE, C( IS, JS ), LDC, B( 1, JS ), LDB, ONE, F( IS, 1 ), LDF )                   CALL SGEMM( 'N', 'T', MB, JS-1, NB, ONE, F( IS, JS ), LDF, E( 1, JS ), LDE, ONE, F( IS, 1 ), LDF )
               }
               if ( I.LT.P ) {
                  CALL SGEMM( 'T', 'N', M-IE, NB, MB, -ONE, A( IS, IE+1 ), LDA, C( IS, JS ), LDC, ONE, C( IE+1, JS ), LDC )                   CALL SGEMM( 'T', 'N', M-IE, NB, MB, -ONE, D( IS, IE+1 ), LDD, F( IS, JS ), LDF, ONE, C( IE+1, JS ), LDC )
               }
  200       CONTINUE
  210    CONTINUE

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of STGSYL

      }
