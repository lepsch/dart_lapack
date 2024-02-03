      SUBROUTINE CCHKBB( NSIZES, MVAL, NVAL, NWDTHS, KK, NTYPES, DOTYPE, NRHS, ISEED, THRESH, NOUNIT, A, LDA, AB, LDAB, BD, BE, Q, LDQ, P, LDP, C, LDC, CC, WORK, LWORK, RWORK, RESULT, INFO )

*  -- LAPACK test routine (input) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDAB, LDC, LDP, LDQ, LWORK, NOUNIT, NRHS, NSIZES, NTYPES, NWDTHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), MVAL( * ), NVAL( * );
      REAL               BD( * ), BE( * ), RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AB( LDAB, * ), C( LDC, * ), CC( LDC, * ), P( LDP, * ), Q( LDQ, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KL, KMAX, KU, M, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CBDT01, CBDT02, CGBBRD, CLACPY, CLASET, CLATMR, CLATMS, CUNT01, SLAHD2, SLASUM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*6, 3*9 /
      DATA               KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3 /
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NTESTT = 0
      INFO = 0

      // Important constants

      BADMM = false;
      BADNN = false;
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      for (J = 1; J <= NSIZES; J++) { // 10
         MMAX = MAX( MMAX, MVAL( J ) )
         IF( MVAL( J ) < 0 ) BADMM = true;
         NMAX = MAX( NMAX, NVAL( J ) )
         IF( NVAL( J ) < 0 ) BADNN = true;
         MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
      } // 10

      BADNNB = false;
      KMAX = 0
      for (J = 1; J <= NWDTHS; J++) { // 20
         KMAX = MAX( KMAX, KK( J ) )
         IF( KK( J ) < 0 ) BADNNB = true;
      } // 20

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1
      } else if ( BADMM ) {
         INFO = -2
      } else if ( BADNN ) {
         INFO = -3
      } else if ( NWDTHS < 0 ) {
         INFO = -4
      } else if ( BADNNB ) {
         INFO = -5
      } else if ( NTYPES < 0 ) {
         INFO = -6
      } else if ( NRHS < 0 ) {
         INFO = -8
      } else if ( LDA < NMAX ) {
         INFO = -13
      } else if ( LDAB < 2*KMAX+1 ) {
         INFO = -15
      } else if ( LDQ < NMAX ) {
         INFO = -19
      } else if ( LDP < NMAX ) {
         INFO = -21
      } else if ( LDC < NMAX ) {
         INFO = -23
      } else if ( ( MAX( LDA, NMAX )+1 )*NMAX.GT.LWORK ) {
         INFO = -26
      }

      if ( INFO != 0 ) {
         xerbla('CCHKBB', -INFO );
         RETURN
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) RETURN;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, widths, types

      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 160
         M = MVAL( JSIZE )
         N = NVAL( JSIZE )
         MNMIN = MIN( M, N )
         AMNINV = ONE / REAL( MAX( 1, M, N ) )

         for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) { // 150
            K = KK( JWIDTH )
            if (K.GE.M && K.GE.N) GO TO 150;
            KL = MAX( 0, MIN( M-1, K ) )
            KU = MAX( 0, MIN( N-1, K ) )

            if ( NSIZES != 1 ) {
               MTYPES = MIN( MAXTYP, NTYPES )
            } else {
               MTYPES = MIN( MAXTYP+1, NTYPES )
            }

            for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 140
               IF( .NOT.DOTYPE( JTYPE ) ) GO TO 140
               NMATS = NMATS + 1
               NTEST = 0

               for (J = 1; J <= 4; J++) { // 30
                  IOLDSD( J ) = ISEED( J )
               } // 30

               // Compute "A".

               // Control parameters:

                   // KMAGN  KMODE        KTYPE
               // =1  O(1)   clustered 1  zero
               // =2  large  clustered 2  identity
               // =3  small  exponential  (none)
               // =4         arithmetic   diagonal, (w/ singular values)
               // =5         random log   (none)
               // =6         random       nonhermitian, w/ singular values
               // =7                      (none)
               // =8                      (none)
               // =9                      random nonhermitian

               if (MTYPES.GT.MAXTYP) GO TO 90;

               ITYPE = KTYPE( JTYPE )
               IMODE = KMODE( JTYPE )

               // Compute norm

               GO TO ( 40, 50, 60 )KMAGN( JTYPE )

               } // 40
               ANORM = ONE
               GO TO 70

               } // 50
               ANORM = ( RTOVFL*ULP )*AMNINV
               GO TO 70

               } // 60
               ANORM = RTUNFL*MAX( M, N )*ULPINV
               GO TO 70

               } // 70

               claset('Full', LDA, N, CZERO, CZERO, A, LDA );
               claset('Full', LDAB, N, CZERO, CZERO, AB, LDAB );
               IINFO = 0
               COND = ULPINV

               // Special Matrices -- Identity & Jordan block

                  // Zero

               if ( ITYPE == 1 ) {
                  IINFO = 0

               } else if ( ITYPE == 2 ) {

                  // Identity

                  for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                     A( JCOL, JCOL ) = ANORM
                  } // 80

               } else if ( ITYPE == 4 ) {

                  // Diagonal Matrix, singular values specified

                  clatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );

               } else if ( ITYPE == 6 ) {

                  // Nonhermitian, singular values specified

                  clatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, KL, KU, 'N', A, LDA, WORK, IINFO );

               } else if ( ITYPE == 9 ) {

                  // Nonhermitian, random entries

                  clatmr(M, N, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, KL, KU, ZERO, ANORM, 'N', A, LDA, IDUMMA, IINFO );

               } else {

                  IINFO = 1
               }

               // Generate Right-Hand Side

               clatmr(M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IDUMMA, M, NRHS, ZERO, ONE, 'NO', C, LDC, IDUMMA, IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               } // 90

               // Copy A to band storage.

               for (J = 1; J <= N; J++) { // 110
                  DO 100 I = MAX( 1, J-KU ), MIN( M, J+KL )
                     AB( KU+1+I-J, J ) = A( I, J )
                  } // 100
               } // 110

               // Copy C

               clacpy('Full', M, NRHS, C, LDC, CC, LDC );

               // Call CGBBRD to compute B, Q and P, and to update C.

               cgbbrd('B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, Q, LDQ, P, LDP, CC, LDC, WORK, RWORK, IINFO );

               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CGBBRD', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO < 0 ) {
                     RETURN
                  } else {
                     RESULT( 1 ) = ULPINV
                     GO TO 120
                  }
               }

               // Test 1:  Check the decomposition A := Q * B * P'
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of P
                    // 4:  Check the computation of Q' * C

               cbdt01(M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RWORK, RESULT( 1 ) );
               cunt01('Columns', M, M, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) );
               cunt01('Rows', N, N, P, LDP, WORK, LWORK, RWORK, RESULT( 3 ) );
               cbdt02(M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RWORK, RESULT( 4 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 4
               } // 120
               NTESTT = NTESTT + NTEST

               // Print out tests which fail.

               for (JR = 1; JR <= NTEST; JR++) { // 130
                  if ( RESULT( JR ).GE.THRESH ) {
                     if (NERRS == 0) CALL SLAHD2( NOUNIT, 'CBB' );
                     NERRS = NERRS + 1
                     WRITE( NOUNIT, FMT = 9998 )M, N, K, IOLDSD, JTYPE, JR, RESULT( JR )
                  }
               } // 130

            } // 140
         } // 150
      } // 160

      // Summary

      slasum('CBB', NOUNIT, NERRS, NTESTT );
      RETURN

 9999 FORMAT( ' CCHKBB: ', A, ' returned INFO=', I5, '.', / 9X, 'M=', I5, ' N=', I5, ' K=', I5, ', JTYPE=', I5, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' M =', I4, ' N=', I4, ', K=', I3, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 )

      // End of CCHKBB

      }
