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
      LQUERY = ( LWORK == -1 )

      if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -1
      } else if ( NOTRAN ) {
         if ( ( IJOB < 0 ) || ( IJOB > 4 ) ) {
            INFO = -2
         }
      }
      if ( INFO == 0 ) {
         if ( M.LE.0 ) {
            INFO = -3
         } else if ( N.LE.0 ) {
            INFO = -4
         } else if ( LDA < MAX( 1, M ) ) {
            INFO = -6
         } else if ( LDB < MAX( 1, N ) ) {
            INFO = -8
         } else if ( LDC < MAX( 1, M ) ) {
            INFO = -10
         } else if ( LDD < MAX( 1, M ) ) {
            INFO = -12
         } else if ( LDE < MAX( 1, N ) ) {
            INFO = -14
         } else if ( LDF < MAX( 1, M ) ) {
            INFO = -16
         }
      }

      if ( INFO == 0 ) {
         if ( NOTRAN ) {
            if ( IJOB == 1 || IJOB == 2 ) {
               LWMIN = MAX( 1, 2*M*N )
            } else {
               LWMIN = 1
            }
         } else {
            LWMIN = 1
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -20
         }
      }

      if ( INFO != 0 ) {
         xerbla('STGSYL', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         SCALE = 1
         if ( NOTRAN ) {
            if ( IJOB != 0 ) {
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
         if ( IJOB >= 3 ) {
            IFUNC = IJOB - 2
            slaset('F', M, N, ZERO, ZERO, C, LDC );
            slaset('F', M, N, ZERO, ZERO, F, LDF );
         } else if ( IJOB >= 1 && NOTRAN ) {
            ISOLVE = 2
         }
      }

      if ( ( MB.LE.1 && NB.LE.1 ) || ( MB >= M && NB >= N ) ) {

         for (IROUND = 1; IROUND <= ISOLVE; IROUND++) { // 30

            // Use unblocked Level 2 solver

            DSCALE = ZERO
            DSUM = ONE
            PQ = 0
            stgsy2(TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, DSUM, DSCALE, IWORK, PQ, INFO );
            if ( DSCALE != ZERO ) {
               if ( IJOB == 1 || IJOB == 3 ) {
                  DIF = SQRT( REAL( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               } else {
                  DIF = SQRT( REAL( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               }
            }

            if ( ISOLVE == 2 && IROUND == 1 ) {
               if ( NOTRAN ) {
                  IFUNC = IJOB
               }
               SCALE2 = SCALE
               slacpy('F', M, N, C, LDC, WORK, M );
               slacpy('F', M, N, F, LDF, WORK( M*N+1 ), M );
               slaset('F', M, N, ZERO, ZERO, C, LDC );
               slaset('F', M, N, ZERO, ZERO, F, LDF );
            } else if ( ISOLVE == 2 && IROUND == 2 ) {
               slacpy('F', M, N, WORK, M, C, LDC );
               slacpy('F', M, N, WORK( M*N+1 ), M, F, LDF );
               SCALE = SCALE2
            }
         } // 30

         RETURN
      }

      // Determine block structure of A

      P = 0
      I = 1
      } // 40
      if (I > M) GO TO 50;
      P = P + 1
      IWORK( P ) = I
      I = I + MB
      if (I >= M) GO TO 50       IF( A( I, I-1 ) != ZERO ) I = I + 1;
      GO TO 40
      } // 50

      IWORK( P+1 ) = M + 1
      IF( IWORK( P ) == IWORK( P+1 ) ) P = P - 1

      // Determine block structure of B

      Q = P + 1
      J = 1
      } // 60
      if (J > N) GO TO 70;
      Q = Q + 1
      IWORK( Q ) = J
      J = J + NB
      if (J >= N) GO TO 70       IF( B( J, J-1 ) != ZERO ) J = J + 1;
      GO TO 60
      } // 70

      IWORK( Q+1 ) = N + 1
      IF( IWORK( Q ) == IWORK( Q+1 ) ) Q = Q - 1

      if ( NOTRAN ) {

         for (IROUND = 1; IROUND <= ISOLVE; IROUND++) { // 150

            // Solve (I, J)-subsystem
                // A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
                // D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
            // for I = P, P - 1,..., 1; J = 1, 2,..., Q

            DSCALE = ZERO
            DSUM = ONE
            PQ = 0
            SCALE = ONE
            for (J = P + 2; J <= Q; J++) { // 130
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               DO 120 I = P, 1, -1
                  IS = IWORK( I )
                  IE = IWORK( I+1 ) - 1
                  MB = IE - IS + 1
                  PPQQ = 0
                  stgsy2(TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, IWORK( Q+2 ), PPQQ, LINFO );
                  if (LINFO > 0) INFO = LINFO;

                  PQ = PQ + PPQQ
                  if ( SCALOC != ONE ) {
                     for (K = 1; K <= JS - 1; K++) { // 80
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 80
                     for (K = JS; K <= JE; K++) { // 90
                        sscal(IS-1, SCALOC, C( 1, K ), 1 );
                        sscal(IS-1, SCALOC, F( 1, K ), 1 );
                     } // 90
                     for (K = JS; K <= JE; K++) { // 100
                        sscal(M-IE, SCALOC, C( IE+1, K ), 1 );
                        sscal(M-IE, SCALOC, F( IE+1, K ), 1 );
                     } // 100
                     for (K = JE + 1; K <= N; K++) { // 110
                        sscal(M, SCALOC, C( 1, K ), 1 );
                        sscal(M, SCALOC, F( 1, K ), 1 );
                     } // 110
                     SCALE = SCALE*SCALOC
                  }

                  // Substitute R(I, J) and L(I, J) into remaining
                  // equation.

                  if ( I > 1 ) {
                     sgemm('N', 'N', IS-1, NB, MB, -ONE, A( 1, IS ), LDA, C( IS, JS ), LDC, ONE, C( 1, JS ), LDC );
                     sgemm('N', 'N', IS-1, NB, MB, -ONE, D( 1, IS ), LDD, C( IS, JS ), LDC, ONE, F( 1, JS ), LDF );
                  }
                  if ( J < Q ) {
                     sgemm('N', 'N', MB, N-JE, NB, ONE, F( IS, JS ), LDF, B( JS, JE+1 ), LDB, ONE, C( IS, JE+1 ), LDC );
                     sgemm('N', 'N', MB, N-JE, NB, ONE, F( IS, JS ), LDF, E( JS, JE+1 ), LDE, ONE, F( IS, JE+1 ), LDF );
                  }
               } // 120
            } // 130
            if ( DSCALE != ZERO ) {
               if ( IJOB == 1 || IJOB == 3 ) {
                  DIF = SQRT( REAL( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               } else {
                  DIF = SQRT( REAL( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               }
            }
            if ( ISOLVE == 2 && IROUND == 1 ) {
               if ( NOTRAN ) {
                  IFUNC = IJOB
               }
               SCALE2 = SCALE
               slacpy('F', M, N, C, LDC, WORK, M );
               slacpy('F', M, N, F, LDF, WORK( M*N+1 ), M );
               slaset('F', M, N, ZERO, ZERO, C, LDC );
               slaset('F', M, N, ZERO, ZERO, F, LDF );
            } else if ( ISOLVE == 2 && IROUND == 2 ) {
               slacpy('F', M, N, WORK, M, C, LDC );
               slacpy('F', M, N, WORK( M*N+1 ), M, F, LDF );
               SCALE = SCALE2
            }
         } // 150

      } else {

         // Solve transposed (I, J)-subsystem
              // A(I, I)**T * R(I, J)  + D(I, I)**T * L(I, J)  =  C(I, J)
              // R(I, J)  * B(J, J)**T + L(I, J)  * E(J, J)**T = -F(I, J)
         // for I = 1,2,..., P; J = Q, Q-1,..., 1

         SCALE = ONE
         for (I = 1; I <= P; I++) { // 210
            IS = IWORK( I )
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 200 J = Q, P + 2, -1
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               stgsy2(TRANS, IFUNC, MB, NB, A( IS, IS ), LDA, B( JS, JS ), LDB, C( IS, JS ), LDC, D( IS, IS ), LDD, E( JS, JS ), LDE, F( IS, JS ), LDF, SCALOC, DSUM, DSCALE, IWORK( Q+2 ), PPQQ, LINFO );
               if (LINFO > 0) INFO = LINFO;
               if ( SCALOC != ONE ) {
                  for (K = 1; K <= JS - 1; K++) { // 160
                     sscal(M, SCALOC, C( 1, K ), 1 );
                     sscal(M, SCALOC, F( 1, K ), 1 );
                  } // 160
                  for (K = JS; K <= JE; K++) { // 170
                     sscal(IS-1, SCALOC, C( 1, K ), 1 );
                     sscal(IS-1, SCALOC, F( 1, K ), 1 );
                  } // 170
                  for (K = JS; K <= JE; K++) { // 180
                     sscal(M-IE, SCALOC, C( IE+1, K ), 1 );
                     sscal(M-IE, SCALOC, F( IE+1, K ), 1 );
                  } // 180
                  for (K = JE + 1; K <= N; K++) { // 190
                     sscal(M, SCALOC, C( 1, K ), 1 );
                     sscal(M, SCALOC, F( 1, K ), 1 );
                  } // 190
                  SCALE = SCALE*SCALOC
               }

               // Substitute R(I, J) and L(I, J) into remaining equation.

               if ( J > P+2 ) {
                  sgemm('N', 'T', MB, JS-1, NB, ONE, C( IS, JS ), LDC, B( 1, JS ), LDB, ONE, F( IS, 1 ), LDF );
                  sgemm('N', 'T', MB, JS-1, NB, ONE, F( IS, JS ), LDF, E( 1, JS ), LDE, ONE, F( IS, 1 ), LDF );
               }
               if ( I < P ) {
                  sgemm('T', 'N', M-IE, NB, MB, -ONE, A( IS, IE+1 ), LDA, C( IS, JS ), LDC, ONE, C( IE+1, JS ), LDC );
                  sgemm('T', 'N', M-IE, NB, MB, -ONE, D( IS, IE+1 ), LDD, F( IS, JS ), LDF, ONE, C( IE+1, JS ), LDC );
               }
            } // 200
         } // 210

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of STGSYL

      }
