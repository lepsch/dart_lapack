      PROGRAM ZBLAT2

*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 17 ;
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      int                NMAX, INCMAX;
      const              NMAX = 65, INCMAX = 2 ;
      int                NINMAX, NIDMAX, NKBMAX, NALMAX, NBEMAX;
      const              NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      double             EPS, ERR, THRESH;
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB, NOUT, NTRA;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR;
      String             TRANS;
      String             SNAMET;
      String             SNAPS, SUMMRY;
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BET( NBEMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( 2*NMAX )
      double             G( NMAX );
      int                IDIM( NIDMAX ), INC( NINMAX ), KB( NKBMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      double             DDIFF;
      bool               LZE;
      // EXTERNAL DDIFF, LZE
      // .. External Subroutines ..
      // EXTERNAL ZCHK1, ZCHK2, ZCHK3, ZCHK4, ZCHK5, ZCHK6, ZCHKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      COMMON             /SRNAMC/SRNAMT
      // .. Data statements ..
      DATA               SNAMES/'ZGEMV ', 'ZGBMV ', 'ZHEMV ', 'ZHBMV ', 'ZHPMV ', 'ZTRMV ', 'ZTBMV ', 'ZTPMV ', 'ZTRSV ', 'ZTBSV ', 'ZTPSV ', 'ZGERC ', 'ZGERU ', 'ZHER  ', 'ZHPR  ', 'ZHER2 ', 'ZHPR2 '/
      // .. Executable Statements ..

      // Read name and unit number for summary output file and open file.

      READ( NIN, FMT = * )SUMMRY
      READ( NIN, FMT = * )NOUT
      OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
      NOUTC = NOUT

      // Read name and unit number for snapshot output file and open file.

      READ( NIN, FMT = * )SNAPS
      READ( NIN, FMT = * )NTRA
      TRACE = NTRA.GE.0
      if ( TRACE ) {
         OPEN( NTRA, FILE = SNAPS, STATUS = 'UNKNOWN' )
      }
      // Read the flag that directs rewinding of the snapshot file.
      READ( NIN, FMT = * )REWI
      REWI = REWI.AND.TRACE
      // Read the flag that directs stopping on any failure.
      READ( NIN, FMT = * )SFATAL
      // Read the flag that indicates whether error exits are to be tested.
      READ( NIN, FMT = * )TSTERR
      // Read the threshold value of the test ratio
      READ( NIN, FMT = * )THRESH

      // Read and check the parameter values for the tests.

      // Values of N
      READ( NIN, FMT = * )NIDIM
      if ( NIDIM.LT.1.OR.NIDIM.GT.NIDMAX ) {
         WRITE( NOUT, FMT = 9997 )'N', NIDMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( IDIM( I ), I = 1, NIDIM )
      DO 10 I = 1, NIDIM
         if ( IDIM( I ).LT.0.OR.IDIM( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX
            GO TO 230
         }
   10 CONTINUE
      // Values of K
      READ( NIN, FMT = * )NKB
      if ( NKB.LT.1.OR.NKB.GT.NKBMAX ) {
         WRITE( NOUT, FMT = 9997 )'K', NKBMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( KB( I ), I = 1, NKB )
      DO 20 I = 1, NKB
         if ( KB( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9995 )
            GO TO 230
         }
   20 CONTINUE
      // Values of INCX and INCY
      READ( NIN, FMT = * )NINC
      if ( NINC.LT.1.OR.NINC.GT.NINMAX ) {
         WRITE( NOUT, FMT = 9997 )'INCX AND INCY', NINMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( INC( I ), I = 1, NINC )
      DO 30 I = 1, NINC
         if ( INC( I ).EQ.0.OR.ABS( INC( I ) ).GT.INCMAX ) {
            WRITE( NOUT, FMT = 9994 )INCMAX
            GO TO 230
         }
   30 CONTINUE
      // Values of ALPHA
      READ( NIN, FMT = * )NALF
      if ( NALF.LT.1.OR.NALF.GT.NALMAX ) {
         WRITE( NOUT, FMT = 9997 )'ALPHA', NALMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( ALF( I ), I = 1, NALF )
      // Values of BETA
      READ( NIN, FMT = * )NBET
      if ( NBET.LT.1.OR.NBET.GT.NBEMAX ) {
         WRITE( NOUT, FMT = 9997 )'BETA', NBEMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( BET( I ), I = 1, NBET )

      // Report values of parameters.

      WRITE( NOUT, FMT = 9993 )
      WRITE( NOUT, FMT = 9992 )( IDIM( I ), I = 1, NIDIM )
      WRITE( NOUT, FMT = 9991 )( KB( I ), I = 1, NKB )
      WRITE( NOUT, FMT = 9990 )( INC( I ), I = 1, NINC )
      WRITE( NOUT, FMT = 9989 )( ALF( I ), I = 1, NALF )
      WRITE( NOUT, FMT = 9988 )( BET( I ), I = 1, NBET )
      if ( .NOT.TSTERR ) {
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9980 )
      }
      WRITE( NOUT, FMT = * )
      WRITE( NOUT, FMT = 9999 )THRESH
      WRITE( NOUT, FMT = * )

      // Read names of subroutines and flags which indicate
      // whether they are to be tested.

      DO 40 I = 1, NSUBS
         LTEST( I ) = .FALSE.
   40 CONTINUE
   50 READ( NIN, FMT = 9984, END = 80 )SNAMET, LTESTT
      DO 60 I = 1, NSUBS
         IF( SNAMET.EQ.SNAMES( I ) ) GO TO 70
   60 CONTINUE
      WRITE( NOUT, FMT = 9986 )SNAMET
      STOP
   70 LTEST( I ) = LTESTT
      GO TO 50

   80 CONTINUE
      CLOSE ( NIN )

      // Compute EPS (the machine precision).

      EPS = EPSILON(RZERO)
      WRITE( NOUT, FMT = 9998 )EPS

      // Check the reliability of ZMVCH using exact data.

      N = MIN( 32, NMAX )
      DO 120 J = 1, N
         DO 110 I = 1, N
            A( I, J ) = MAX( I - J + 1, 0 )
  110    CONTINUE
         X( J ) = J
         Y( J ) = ZERO
  120 CONTINUE
      DO 130 J = 1, N
         YY( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
  130 CONTINUE
      // YY holds the exact result. On exit from ZMVCH YT holds
      // the result computed by ZMVCH.
      TRANS = 'N'
      zmvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LZE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }
      TRANS = 'T'
      zmvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
      SAME = LZE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR.NE.RZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }

      // Test each subroutine in turn.

      DO 210 ISNUM = 1, NSUBS
         WRITE( NOUT, FMT = * )
         if ( .NOT.LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM )
         } else {
            SRNAMT = SNAMES( ISNUM )
            // Test error exits.
            if ( TSTERR ) {
               zchke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * )
            }
            // Test computations.
            INFOT = 0
            OK = .TRUE.
            FATAL = .FALSE.
            GO TO ( 140, 140, 150, 150, 150, 160, 160, 160, 160, 160, 160, 170, 170, 180, 180, 190, 190 )ISNUM
            // Test ZGEMV, 01, and ZGBMV, 02.
  140       CALL ZCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
            // Test ZHEMV, 03, ZHBMV, 04, and ZHPMV, 05.
  150       CALL ZCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
            // Test ZTRMV, 06, ZTBMV, 07, ZTPMV, 08,
            // ZTRSV, 09, ZTBSV, 10, and ZTPSV, 11.
  160       CALL ZCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test ZGERC, 12, ZGERU, 13.
  170       CALL ZCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test ZHER, 14, and ZHPR, 15.
  180       CALL ZCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test ZHER2, 16, and ZHPR2, 17.
  190       CALL ZCHK6( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

  200       IF( FATAL.AND.SFATAL ) GO TO 220
         }
  210 CONTINUE
      WRITE( NOUT, FMT = 9982 )
      GO TO 240

  220 CONTINUE
      WRITE( NOUT, FMT = 9981 )
      GO TO 240

  230 CONTINUE
      WRITE( NOUT, FMT = 9987 )

  240 CONTINUE
      IF( TRACE ) CLOSE ( NTRA )
      CLOSE ( NOUT )
      STOP

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 )
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, D9.1 )
 9997 FORMAT( ' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 )
 9995 FORMAT( ' VALUE OF K IS LESS THAN 0' )
 9994 FORMAT( ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ', I2 )
 9993 FORMAT( ' TESTS OF THE COMPLEX*16       LEVEL 2 BLAS', //' THE F', 'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9992 FORMAT( '   FOR N              ', 9I6 )
 9991 FORMAT( '   FOR K              ', 7I6 )
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 )
 9989 FORMAT( '   FOR ALPHA          ', 7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9988 FORMAT( '   FOR BETA           ', 7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9987 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' )
 9986 FORMAT( ' SUBPROGRAM NAME ', A6, ' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' )
 9985 FORMAT( ' ERROR IN ZMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' ZMVCH WAS CALLED WITH TRANS = ', A1, ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '.', / ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.' , /' ******* TESTS ABANDONED *******' )
 9984 FORMAT( A6, L2 )
 9983 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9982 FORMAT( /' END OF TESTS' )
 9981 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' )
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )

      // End of ZBLAT2

      }
      SUBROUTINE ZCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )

*  Tests ZGEMV and ZGBMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, BETA, BLS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA, LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK, NL, NS;
      bool               BANDED, FULL, NULL, RESET, SAME, TRAN;
      String             TRANS, TRANSS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZGBMV, ZGEMV, ZMAKE, ZMVCH, ZREGR1
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'NTC'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ).EQ.'E'
      BANDED = SNAME( 3: 3 ).EQ.'B'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 11
      } else if ( BANDED ) {
         NARGS = 13
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 120 IN = 1, NIDIM
         N = IDIM( IN )
         ND = N/2 + 1

         DO 110 IM = 1, 2
            IF( IM.EQ.1 ) M = MAX( N - ND, 0 )             IF( IM.EQ.2 ) M = MIN( N + ND, NMAX )

            if ( BANDED ) {
               NK = NKB
            } else {
               NK = 1
            }
            DO 100 IKU = 1, NK
               if ( BANDED ) {
                  KU = KB( IKU )
                  KL = MAX( KU - 1, 0 )
               } else {
                  KU = N - 1
                  KL = M - 1
               }
               // Set LDA to 1 more than minimum value if room.
               if ( BANDED ) {
                  LDA = KL + KU + 1
               } else {
                  LDA = M
               }
               IF( LDA.LT.NMAX ) LDA = LDA + 1
               // Skip tests if not enough room.
               IF( LDA.GT.NMAX ) GO TO 100
               LAA = LDA*N
               NULL = N.LE.0.OR.M.LE.0

               // Generate the matrix A.

               TRANSL = ZERO
               zmake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL );

               DO 90 IC = 1, 3
                  TRANS = ICH( IC: IC )
                  TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'

                  if ( TRAN ) {
                     ML = N
                     NL = M
                  } else {
                     ML = M
                     NL = N
                  }

                  DO 80 IX = 1, NINC
                     INCX = INC( IX )
                     LX = ABS( INCX )*NL

                     // Generate the vector X.

                     TRANSL = HALF
                     zmake('GE', ' ', ' ', 1, NL, X, 1, XX, ABS( INCX ), 0, NL - 1, RESET, TRANSL );
                     if ( NL.GT.1 ) {
                        X( NL/2 ) = ZERO
                        XX( 1 + ABS( INCX )*( NL/2 - 1 ) ) = ZERO
                     }

                     DO 70 IY = 1, NINC
                        INCY = INC( IY )
                        LY = ABS( INCY )*ML

                        DO 60 IA = 1, NALF
                           ALPHA = ALF( IA )

                           DO 50 IB = 1, NBET
                              BETA = BET( IB )

                              // Generate the vector Y.

                              TRANSL = ZERO
                              zmake('GE', ' ', ' ', 1, ML, Y, 1, YY, ABS( INCY ), 0, ML - 1, RESET, TRANSL );

                              NC = NC + 1

                              // Save every datum before calling the
                              // subroutine.

                              TRANSS = TRANS
                              MS = M
                              NS = N
                              KLS = KL
                              KUS = KU
                              ALS = ALPHA
                              DO 10 I = 1, LAA
                                 AS( I ) = AA( I )
   10                         CONTINUE
                              LDAS = LDA
                              DO 20 I = 1, LX
                                 XS( I ) = XX( I )
   20                         CONTINUE
                              INCXS = INCX
                              BLS = BETA
                              DO 30 I = 1, LY
                                 YS( I ) = YY( I )
   30                         CONTINUE
                              INCYS = INCY

                              // Call the subroutine.

                              if ( FULL ) {
                                 IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
                                 IF( REWI ) REWIND NTRA                                  CALL ZGEMV( TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                              } else if ( BANDED ) {
                                 IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
                                 IF( REWI ) REWIND NTRA                                  CALL ZGBMV( TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                              }

                              // Check if error-exit was taken incorrectly.

                              if ( .NOT.OK ) {
                                 WRITE( NOUT, FMT = 9993 )
                                 FATAL = .TRUE.
                                 GO TO 130
                              }

                              // See what data changed inside subroutines.

                              ISAME( 1 ) = TRANS.EQ.TRANSS
                              ISAME( 2 ) = MS.EQ.M
                              ISAME( 3 ) = NS.EQ.N
                              if ( FULL ) {
                                 ISAME( 4 ) = ALS.EQ.ALPHA
                                 ISAME( 5 ) = LZE( AS, AA, LAA )
                                 ISAME( 6 ) = LDAS.EQ.LDA
                                 ISAME( 7 ) = LZE( XS, XX, LX )
                                 ISAME( 8 ) = INCXS.EQ.INCX
                                 ISAME( 9 ) = BLS.EQ.BETA
                                 if ( NULL ) {
                                    ISAME( 10 ) = LZE( YS, YY, LY )
                                 } else {
                                    ISAME( 10 ) = LZERES( 'GE', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 11 ) = INCYS.EQ.INCY
                              } else if ( BANDED ) {
                                 ISAME( 4 ) = KLS.EQ.KL
                                 ISAME( 5 ) = KUS.EQ.KU
                                 ISAME( 6 ) = ALS.EQ.ALPHA
                                 ISAME( 7 ) = LZE( AS, AA, LAA )
                                 ISAME( 8 ) = LDAS.EQ.LDA
                                 ISAME( 9 ) = LZE( XS, XX, LX )
                                 ISAME( 10 ) = INCXS.EQ.INCX
                                 ISAME( 11 ) = BLS.EQ.BETA
                                 if ( NULL ) {
                                    ISAME( 12 ) = LZE( YS, YY, LY )
                                 } else {
                                    ISAME( 12 ) = LZERES( 'GE', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 13 ) = INCYS.EQ.INCY
                              }

                              // If data was incorrectly changed, report
                              // and return.

                              SAME = .TRUE.
                              DO 40 I = 1, NARGS
                                 SAME = SAME.AND.ISAME( I )
                                 IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                         CONTINUE
                              if ( .NOT.SAME ) {
                                 FATAL = .TRUE.
                                 GO TO 130
                              }

                              if ( .NOT.NULL ) {

                                 // Check the result.

                                 zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
                                 ERRMAX = MAX( ERRMAX, ERR )
                                 // If got really bad answer, report and
                                 // return.
                                 IF( FATAL ) GO TO 130
                              } else {
                                 // Avoid repeating tests with M.le.0 or
                                 // N.le.0.
                                 GO TO 110
                              }

   50                      CONTINUE

   60                   CONTINUE

   70                CONTINUE

   80             CONTINUE

   90          CONTINUE

  100       CONTINUE

  110    CONTINUE

  120 CONTINUE

      // Regression test to verify preservation of y when m zero, n nonzero.

      zregr1(TRANS, M, N, LY, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY, YS );
      if ( FULL ) {
         IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
         IF( REWI ) REWIND NTRA          CALL ZGEMV( TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
      } else if ( BANDED ) {
         IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
         IF( REWI ) REWIND NTRA          CALL ZGBMV( TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
      }
      NC = NC + 1
      if ( .NOT.LZE( YS, YY, LY ) ) {
         WRITE( NOUT, FMT = 9998 )NARGS - 1
         FATAL = .TRUE.
         GO TO 130
      }

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 140

  130 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
      }

  140 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 4( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ') .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ')         .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK1

      }
      SUBROUTINE ZCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )

*  Tests ZHEMV, ZHBMV and ZHPMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, BETA, BLS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IB, IC, IK, IN, INCX, INCXS, INCY, INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHBMV, ZHEMV, ZHPMV, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ).EQ.'E'
      BANDED = SNAME( 3: 3 ).EQ.'B'
      PACKED = SNAME( 3: 3 ).EQ.'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 10
      } else if ( BANDED ) {
         NARGS = 11
      } else if ( PACKED ) {
         NARGS = 9
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 110 IN = 1, NIDIM
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         DO 100 IK = 1, NK
            if ( BANDED ) {
               K = KB( IK )
            } else {
               K = N - 1
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1
            } else {
               LDA = N
            }
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 100
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            DO 90 IC = 1, 2
               UPLO = ICH( IC: IC )

               // Generate the matrix A.

               TRANSL = ZERO
               zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

               DO 80 IX = 1, NINC
                  INCX = INC( IX )
                  LX = ABS( INCX )*N

                  // Generate the vector X.

                  TRANSL = HALF
                  zmake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     X( N/2 ) = ZERO
                     XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                  }

                  DO 70 IY = 1, NINC
                     INCY = INC( IY )
                     LY = ABS( INCY )*N

                     DO 60 IA = 1, NALF
                        ALPHA = ALF( IA )

                        DO 50 IB = 1, NBET
                           BETA = BET( IB )

                           // Generate the vector Y.

                           TRANSL = ZERO
                           zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           UPLOS = UPLO
                           NS = N
                           KS = K
                           ALS = ALPHA
                           DO 10 I = 1, LAA
                              AS( I ) = AA( I )
   10                      CONTINUE
                           LDAS = LDA
                           DO 20 I = 1, LX
                              XS( I ) = XX( I )
   20                      CONTINUE
                           INCXS = INCX
                           BLS = BETA
                           DO 30 I = 1, LY
                              YS( I ) = YY( I )
   30                      CONTINUE
                           INCYS = INCY

                           // Call the subroutine.

                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY
                              IF( REWI ) REWIND NTRA                               CALL ZHEMV( UPLO, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY
                              IF( REWI ) REWIND NTRA                               CALL ZHBMV( UPLO, N, K, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY
                              IF( REWI ) REWIND NTRA                               CALL ZHPMV( UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY )
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( .NOT.OK ) {
                              WRITE( NOUT, FMT = 9992 )
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = UPLO.EQ.UPLOS
                           ISAME( 2 ) = NS.EQ.N
                           if ( FULL ) {
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LZE( AS, AA, LAA )
                              ISAME( 5 ) = LDAS.EQ.LDA
                              ISAME( 6 ) = LZE( XS, XX, LX )
                              ISAME( 7 ) = INCXS.EQ.INCX
                              ISAME( 8 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 9 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 9 ) = LZERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 10 ) = INCYS.EQ.INCY
                           } else if ( BANDED ) {
                              ISAME( 3 ) = KS.EQ.K
                              ISAME( 4 ) = ALS.EQ.ALPHA
                              ISAME( 5 ) = LZE( AS, AA, LAA )
                              ISAME( 6 ) = LDAS.EQ.LDA
                              ISAME( 7 ) = LZE( XS, XX, LX )
                              ISAME( 8 ) = INCXS.EQ.INCX
                              ISAME( 9 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 10 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 10 ) = LZERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 11 ) = INCYS.EQ.INCY
                           } else if ( PACKED ) {
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LZE( AS, AA, LAA )
                              ISAME( 5 ) = LZE( XS, XX, LX )
                              ISAME( 6 ) = INCXS.EQ.INCX
                              ISAME( 7 ) = BLS.EQ.BETA
                              if ( NULL ) {
                                 ISAME( 8 ) = LZE( YS, YY, LY )
                              } else {
                                 ISAME( 8 ) = LZERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 9 ) = INCYS.EQ.INCY
                           }

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = .TRUE.
                           DO 40 I = 1, NARGS
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                      CONTINUE
                           if ( .NOT.SAME ) {
                              FATAL = .TRUE.
                              GO TO 120
                           }

                           if ( .NOT.NULL ) {

                              // Check the result.

                              zmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, .TRUE. );
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              IF( FATAL ) GO TO 120
                           } else {
                              // Avoid repeating tests with N.le.0
                              GO TO 110
                           }

   50                   CONTINUE

   60                CONTINUE

   70             CONTINUE

   80          CONTINUE

   90       CONTINUE

  100    CONTINUE

  110 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

  120 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY
      }

  130 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',(', F4.1, ',', F4.1, '), AP, X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ')                .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), A,', I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2, ')         .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',(', F4.1, ',', F4.1, '), A,', I3, ', X,', I2, ',(', F4.1, ',', F4.1, '), ', 'Y,', I2, ')             .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK2

      }
      SUBROUTINE ZCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, XT, G, Z )

*  Tests ZTRMV, ZTBMV, ZTPMV, ZTRSV, ZTBSV and ZTPSV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XT( NMAX ), XX( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      COMPLEX*16         TRANSL
      double             ERR, ERRMAX;
      int                I, ICD, ICT, ICU, IK, IN, INCX, INCXS, IX, K, KS, LAA, LDA, LDAS, LX, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             DIAG, DIAGS, TRANS, TRANSS, UPLO, UPLOS;
      String             ICHD, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZMAKE, ZMVCH, ZTBMV, ZTBSV, ZTPMV, ZTPSV, ZTRMV, ZTRSV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ).EQ.'R'
      BANDED = SNAME( 3: 3 ).EQ.'B'
      PACKED = SNAME( 3: 3 ).EQ.'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 8
      } else if ( BANDED ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 7
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
      // Set up zero vector for ZMVCH.
      DO 10 I = 1, NMAX
         Z( I ) = ZERO
   10 CONTINUE

      DO 110 IN = 1, NIDIM
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         DO 100 IK = 1, NK
            if ( BANDED ) {
               K = KB( IK )
            } else {
               K = N - 1
            }
            // Set LDA to 1 more than minimum value if room.
            if ( BANDED ) {
               LDA = K + 1
            } else {
               LDA = N
            }
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 100
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            DO 90 ICU = 1, 2
               UPLO = ICHU( ICU: ICU )

               DO 80 ICT = 1, 3
                  TRANS = ICHT( ICT: ICT )

                  DO 70 ICD = 1, 2
                     DIAG = ICHD( ICD: ICD )

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 2: 3 ), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

                     DO 60 IX = 1, NINC
                        INCX = INC( IX )
                        LX = ABS( INCX )*N

                        // Generate the vector X.

                        TRANSL = HALF
                        zmake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
                        if ( N.GT.1 ) {
                           X( N/2 ) = ZERO
                           XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                        }

                        NC = NC + 1

                        // Save every datum before calling the subroutine.

                        UPLOS = UPLO
                        TRANSS = TRANS
                        DIAGS = DIAG
                        NS = N
                        KS = K
                        DO 20 I = 1, LAA
                           AS( I ) = AA( I )
   20                   CONTINUE
                        LDAS = LDA
                        DO 30 I = 1, LX
                           XS( I ) = XX( I )
   30                   CONTINUE
                        INCXS = INCX

                        // Call the subroutine.

                        if ( SNAME( 4: 5 ).EQ.'MV' ) {
                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTRMV( UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTBMV( UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTPMV( UPLO, TRANS, DIAG, N, AA, XX, INCX )
                           }
                        } else if ( SNAME( 4: 5 ).EQ.'SV' ) {
                           if ( FULL ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTRSV( UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX )
                           } else if ( BANDED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTBSV( UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX )
                           } else if ( PACKED ) {
                              IF( TRACE ) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX
                              IF( REWI ) REWIND NTRA                               CALL ZTPSV( UPLO, TRANS, DIAG, N, AA, XX, INCX )
                           }
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9992 )
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLO.EQ.UPLOS
                        ISAME( 2 ) = TRANS.EQ.TRANSS
                        ISAME( 3 ) = DIAG.EQ.DIAGS
                        ISAME( 4 ) = NS.EQ.N
                        if ( FULL ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA )
                           ISAME( 6 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 7 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 7 ) = LZERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 8 ) = INCXS.EQ.INCX
                        } else if ( BANDED ) {
                           ISAME( 5 ) = KS.EQ.K
                           ISAME( 6 ) = LZE( AS, AA, LAA )
                           ISAME( 7 ) = LDAS.EQ.LDA
                           if ( NULL ) {
                              ISAME( 8 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 8 ) = LZERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 9 ) = INCXS.EQ.INCX
                        } else if ( PACKED ) {
                           ISAME( 5 ) = LZE( AS, AA, LAA )
                           if ( NULL ) {
                              ISAME( 6 ) = LZE( XS, XX, LX )
                           } else {
                              ISAME( 6 ) = LZERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 7 ) = INCXS.EQ.INCX
                        }

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = .TRUE.
                        DO 40 I = 1, NARGS
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                   CONTINUE
                        if ( .NOT.SAME ) {
                           FATAL = .TRUE.
                           GO TO 120
                        }

                        if ( .NOT.NULL ) {
                           if ( SNAME( 4: 5 ).EQ.'MV' ) {

                              // Check the result.

                              zmvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, .TRUE. );
                           } else if ( SNAME( 4: 5 ).EQ.'SV' ) {

                              // Compute approximation to original vector.

                              DO 50 I = 1, N
                                 Z( I ) = XX( 1 + ( I - 1 )* ABS( INCX ) )                                  XX( 1 + ( I - 1 )*ABS( INCX ) ) = X( I )
   50                         CONTINUE
                              zmvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, .FALSE. );
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 120
                        } else {
                           // Avoid repeating tests with N.le.0.
                           GO TO 110
                        }

   60                CONTINUE

   70             CONTINUE

   80          CONTINUE

   90       CONTINUE

  100    CONTINUE

  110 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

  120 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX
      }

  130 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', AP, ', 'X,', I2, ')                                      .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), 2( I3, ',' ), ' A,', I3, ', X,', I2, ')                               .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', A,', I3, ', X,', I2, ')                                   .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK3

      }
      SUBROUTINE ZCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests ZGERC and ZGERU.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS, NC, ND, NS;
      bool               CONJ, NULL, RESET, SAME;
      // .. Local Arrays ..
      COMPLEX*16         W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZGERC, ZGERU, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      CONJ = SNAME( 5: 5 ).EQ.'C'
      // Define the number of arguments.
      NARGS = 9

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 120 IN = 1, NIDIM
         N = IDIM( IN )
         ND = N/2 + 1

         DO 110 IM = 1, 2
            IF( IM.EQ.1 ) M = MAX( N - ND, 0 )             IF( IM.EQ.2 ) M = MIN( N + ND, NMAX )

            // Set LDA to 1 more than minimum value if room.
            LDA = M
            IF( LDA.LT.NMAX ) LDA = LDA + 1
            // Skip tests if not enough room.
            IF( LDA.GT.NMAX ) GO TO 110
            LAA = LDA*N
            NULL = N.LE.0.OR.M.LE.0

            DO 100 IX = 1, NINC
               INCX = INC( IX )
               LX = ABS( INCX )*M

               // Generate the vector X.

               TRANSL = HALF
               zmake('GE', ' ', ' ', 1, M, X, 1, XX, ABS( INCX ), 0, M - 1, RESET, TRANSL );
               if ( M.GT.1 ) {
                  X( M/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( M/2 - 1 ) ) = ZERO
               }

               DO 90 IY = 1, NINC
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  DO 80 IA = 1, NALF
                     ALPHA = ALF( IA )

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, M - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     MS = M
                     NS = N
                     ALS = ALPHA
                     DO 10 I = 1, LAA
                        AS( I ) = AA( I )
   10                CONTINUE
                     LDAS = LDA
                     DO 20 I = 1, LX
                        XS( I ) = XX( I )
   20                CONTINUE
                     INCXS = INCX
                     DO 30 I = 1, LY
                        YS( I ) = YY( I )
   30                CONTINUE
                     INCYS = INCY

                     // Call the subroutine.

                     IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA
                     if ( CONJ ) {
                        IF( REWI ) REWIND NTRA                         CALL ZGERC( M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     } else {
                        IF( REWI ) REWIND NTRA                         CALL ZGERU( M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9993 )
                        FATAL = .TRUE.
                        GO TO 140
                     }

                     // See what data changed inside subroutine.

                     ISAME( 1 ) = MS.EQ.M
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LZE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LZE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LZERES( 'GE', ' ', M, N, AS, AA, LDA )
                     }
                     ISAME( 9 ) = LDAS.EQ.LDA

                     // If data was incorrectly changed, report and return.

                     SAME = .TRUE.
                     DO 40 I = 1, NARGS
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                CONTINUE
                     if ( .NOT.SAME ) {
                        FATAL = .TRUE.
                        GO TO 140
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           DO 50 I = 1, M
                              Z( I ) = X( I )
   50                      CONTINUE
                        } else {
                           DO 60 I = 1, M
                              Z( I ) = X( M - I + 1 )
   60                      CONTINUE
                        }
                        DO 70 J = 1, N
                           if ( INCY.GT.0 ) {
                              W( 1 ) = Y( J )
                           } else {
                              W( 1 ) = Y( N - J + 1 )
                           }
                           IF( CONJ ) W( 1 ) = DCONJG( W( 1 ) )                            CALL ZMVCH( 'N', M, 1, ALPHA, Z, NMAX, W, 1, ONE, A( 1, J ), 1, YT, G, AA( 1 + ( J - 1 )*LDA ), EPS, ERR, FATAL, NOUT, .TRUE. )
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 130
   70                   CONTINUE
                     } else {
                        // Avoid repeating tests with M.le.0 or N.le.0.
                        GO TO 110
                     }

   80             CONTINUE

   90          CONTINUE

  100       CONTINUE

  110    CONTINUE

  120 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 150

  130 CONTINUE
      WRITE( NOUT, FMT = 9995 )J

  140 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA

  150 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( I3, ',' ), '(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ')                   ', '      .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK4

      }
      SUBROUTINE ZCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests ZHER and ZHPR.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, TRANSL
      double             ERR, ERRMAX, RALPHA, RALS;
      int                I, IA, IC, IN, INCX, INCXS, IX, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      COMPLEX*16         W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHER, ZHPR, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ).EQ.'E'
      PACKED = SNAME( 3: 3 ).EQ.'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 7
      } else if ( PACKED ) {
         NARGS = 6
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 100 IN = 1, NIDIM
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         IF( LDA.LT.NMAX ) LDA = LDA + 1
         // Skip tests if not enough room.
         IF( LDA.GT.NMAX ) GO TO 100
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         DO 90 IC = 1, 2
            UPLO = ICH( IC: IC )
            UPPER = UPLO.EQ.'U'

            DO 80 IX = 1, NINC
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               zmake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               DO 70 IA = 1, NALF
                  RALPHA = DBLE( ALF( IA ) )
                  ALPHA = DCMPLX( RALPHA, RZERO )
                  NULL = N.LE.0.OR.RALPHA.EQ.RZERO

                  // Generate the matrix A.

                  TRANSL = ZERO
                  zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                  NC = NC + 1

                  // Save every datum before calling the subroutine.

                  UPLOS = UPLO
                  NS = N
                  RALS = RALPHA
                  DO 10 I = 1, LAA
                     AS( I ) = AA( I )
   10             CONTINUE
                  LDAS = LDA
                  DO 20 I = 1, LX
                     XS( I ) = XX( I )
   20             CONTINUE
                  INCXS = INCX

                  // Call the subroutine.

                  if ( FULL ) {
                     IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, RALPHA, INCX, LDA
                     IF( REWI ) REWIND NTRA
                     zher(UPLO, N, RALPHA, XX, INCX, AA, LDA );
                  } else if ( PACKED ) {
                     IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, RALPHA, INCX
                     IF( REWI ) REWIND NTRA
                     zhpr(UPLO, N, RALPHA, XX, INCX, AA );
                  }

                  // Check if error-exit was taken incorrectly.

                  if ( .NOT.OK ) {
                     WRITE( NOUT, FMT = 9992 )
                     FATAL = .TRUE.
                     GO TO 120
                  }

                  // See what data changed inside subroutines.

                  ISAME( 1 ) = UPLO.EQ.UPLOS
                  ISAME( 2 ) = NS.EQ.N
                  ISAME( 3 ) = RALS.EQ.RALPHA
                  ISAME( 4 ) = LZE( XS, XX, LX )
                  ISAME( 5 ) = INCXS.EQ.INCX
                  if ( NULL ) {
                     ISAME( 6 ) = LZE( AS, AA, LAA )
                  } else {
                     ISAME( 6 ) = LZERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA )
                  }
                  if ( .NOT.PACKED ) {
                     ISAME( 7 ) = LDAS.EQ.LDA
                  }

                  // If data was incorrectly changed, report and return.

                  SAME = .TRUE.
                  DO 30 I = 1, NARGS
                     SAME = SAME.AND.ISAME( I )
                     IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   30             CONTINUE
                  if ( .NOT.SAME ) {
                     FATAL = .TRUE.
                     GO TO 120
                  }

                  if ( .NOT.NULL ) {

                     // Check the result column by column.

                     if ( INCX.GT.0 ) {
                        DO 40 I = 1, N
                           Z( I ) = X( I )
   40                   CONTINUE
                     } else {
                        DO 50 I = 1, N
                           Z( I ) = X( N - I + 1 )
   50                   CONTINUE
                     }
                     JA = 1
                     DO 60 J = 1, N
                        W( 1 ) = DCONJG( Z( J ) )
                        if ( UPPER ) {
                           JJ = 1
                           LJ = J
                        } else {
                           JJ = J
                           LJ = N - J + 1
                        }
                        zmvch('N', LJ, 1, ALPHA, Z( JJ ), LJ, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, .TRUE. );
                        if ( FULL ) {
                           if ( UPPER ) {
                              JA = JA + LDA
                           } else {
                              JA = JA + LDA + 1
                           }
                        } else {
                           JA = JA + LJ
                        }
                        ERRMAX = MAX( ERRMAX, ERR )
                        // If got really bad answer, report and return.
                        IF( FATAL ) GO TO 110
   60                CONTINUE
                  } else {
                     // Avoid repeating tests if N.le.0.
                     IF( N.LE.0 ) GO TO 100
                  }

   70          CONTINUE

   80       CONTINUE

   90    CONTINUE

  100 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

  110 CONTINUE
      WRITE( NOUT, FMT = 9995 )J

  120 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, RALPHA, INCX, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, RALPHA, INCX
      }

  130 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', AP)                                         .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', A,', I3, ')                                      .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK5

      }
      SUBROUTINE ZCHK6( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests ZHER2 and ZHPR2.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, HALF, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), HALF = ( 0.5D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      // .. Scalar Arguments ..
      double             EPS, THRESH;
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX, 2 )
      double             G( NMAX );
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, ALS, TRANSL
      double             ERR, ERRMAX;
      int                I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      COMPLEX*16         W( 2 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LZE, LZERES;
      // EXTERNAL LZE, LZERES
      // .. External Subroutines ..
      // EXTERNAL ZHER2, ZHPR2, ZMAKE, ZMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ).EQ.'E'
      PACKED = SNAME( 3: 3 ).EQ.'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 8
      }

      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO

      DO 140 IN = 1, NIDIM
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         IF( LDA.LT.NMAX ) LDA = LDA + 1
         // Skip tests if not enough room.
         IF( LDA.GT.NMAX ) GO TO 140
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         DO 130 IC = 1, 2
            UPLO = ICH( IC: IC )
            UPPER = UPLO.EQ.'U'

            DO 120 IX = 1, NINC
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               zmake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               DO 110 IY = 1, NINC
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  zmake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  DO 100 IA = 1, NALF
                     ALPHA = ALF( IA )
                     NULL = N.LE.0.OR.ALPHA.EQ.ZERO

                     // Generate the matrix A.

                     TRANSL = ZERO
                     zmake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     UPLOS = UPLO
                     NS = N
                     ALS = ALPHA
                     DO 10 I = 1, LAA
                        AS( I ) = AA( I )
   10                CONTINUE
                     LDAS = LDA
                     DO 20 I = 1, LX
                        XS( I ) = XX( I )
   20                CONTINUE
                     INCXS = INCX
                     DO 30 I = 1, LY
                        YS( I ) = YY( I )
   30                CONTINUE
                     INCYS = INCY

                     // Call the subroutine.

                     if ( FULL ) {
                        IF( TRACE ) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA
                        IF( REWI ) REWIND NTRA                         CALL ZHER2( UPLO, N, ALPHA, XX, INCX, YY, INCY, AA, LDA )
                     } else if ( PACKED ) {
                        IF( TRACE ) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY
                        IF( REWI ) REWIND NTRA                         CALL ZHPR2( UPLO, N, ALPHA, XX, INCX, YY, INCY, AA )
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9992 )
                        FATAL = .TRUE.
                        GO TO 160
                     }

                     // See what data changed inside subroutines.

                     ISAME( 1 ) = UPLO.EQ.UPLOS
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LZE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LZE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LZE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LZERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA )
                     }
                     if ( .NOT.PACKED ) {
                        ISAME( 9 ) = LDAS.EQ.LDA
                     }

                     // If data was incorrectly changed, report and return.

                     SAME = .TRUE.
                     DO 40 I = 1, NARGS
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
   40                CONTINUE
                     if ( .NOT.SAME ) {
                        FATAL = .TRUE.
                        GO TO 160
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           DO 50 I = 1, N
                              Z( I, 1 ) = X( I )
   50                      CONTINUE
                        } else {
                           DO 60 I = 1, N
                              Z( I, 1 ) = X( N - I + 1 )
   60                      CONTINUE
                        }
                        if ( INCY.GT.0 ) {
                           DO 70 I = 1, N
                              Z( I, 2 ) = Y( I )
   70                      CONTINUE
                        } else {
                           DO 80 I = 1, N
                              Z( I, 2 ) = Y( N - I + 1 )
   80                      CONTINUE
                        }
                        JA = 1
                        DO 90 J = 1, N
                           W( 1 ) = ALPHA*DCONJG( Z( J, 2 ) )
                           W( 2 ) = DCONJG( ALPHA )*DCONJG( Z( J, 1 ) )
                           if ( UPPER ) {
                              JJ = 1
                              LJ = J
                           } else {
                              JJ = J
                              LJ = N - J + 1
                           }
                           zmvch('N', LJ, 2, ONE, Z( JJ, 1 ), NMAX, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, .TRUE. );
                           if ( FULL ) {
                              if ( UPPER ) {
                                 JA = JA + LDA
                              } else {
                                 JA = JA + LDA + 1
                              }
                           } else {
                              JA = JA + LJ
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           IF( FATAL ) GO TO 150
   90                   CONTINUE
                     } else {
                        // Avoid repeating tests with N.le.0.
                        IF( N.LE.0 ) GO TO 140
                     }

  100             CONTINUE

  110          CONTINUE

  120       CONTINUE

  130    CONTINUE

  140 CONTINUE

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 170

  150 CONTINUE
      WRITE( NOUT, FMT = 9995 )J

  160 CONTINUE
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY
      }

  170 CONTINUE
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', AP)                     ', '       .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',(', F4.1, ',', F4.1, '), X,', I2, ', Y,', I2, ', A,', I3, ')             ', '            .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of ZCHK6

      }
      SUBROUTINE ZCHKE( ISNUM, SRNAMT, NOUT )

*  Tests the error exits from the Level 2 Blas.
*  Requires a special version of the error-handling routine XERBLA.
*  ALPHA, RALPHA, BETA, A, X and Y should not need to be defined.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                ISNUM, NOUT;
      String             SRNAMT;
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, BETA
      double             RALPHA;
      // .. Local Arrays ..
      COMPLEX*16         A( 1, 1 ), X( 1 ), Y( 1 )
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGBMV, ZGEMV, ZGERC, ZGERU, ZHBMV, ZHEMV, ZHER, ZHER2, ZHPMV, ZHPR, ZHPR2, ZTBMV, ZTBSV, ZTPMV, ZTPSV, ZTRMV, ZTRSV
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      // OK is set to .FALSE. by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      OK = .TRUE.
      // LERR is set to .TRUE. by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      LERR = .FALSE.
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170 )ISNUM
   10 INFOT = 1
      zgemv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgemv('N', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      zgemv('N', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      zgemv('N', 2, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      zgemv('N', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      zgemv('N', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   20 INFOT = 1
      zgbmv('/', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgbmv('N', -1, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      zgbmv('N', 0, -1, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgbmv('N', 0, 0, -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgbmv('N', 2, 0, 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      zgbmv('N', 0, 0, 1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   30 INFOT = 1
      zhemv('/', 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhemv('U', -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zhemv('U', 2, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zhemv('U', 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      zhemv('U', 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   40 INFOT = 1
      zhbmv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhbmv('U', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      zhbmv('U', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      zhbmv('U', 0, 1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   50 INFOT = 1
      zhpmv('/', 0, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhpmv('U', -1, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      zhpmv('U', 0, ALPHA, A, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      zhpmv('U', 0, ALPHA, A, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   60 INFOT = 1
      ztrmv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztrmv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztrmv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztrmv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      ztrmv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      ztrmv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   70 INFOT = 1
      ztbmv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztbmv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztbmv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztbmv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      ztbmv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ztbmv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      ztbmv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   80 INFOT = 1
      ztpmv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztpmv('U', '/', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztpmv('U', 'N', '/', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztpmv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ztpmv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
   90 INFOT = 1
      ztrsv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztrsv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztrsv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztrsv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      ztrsv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      ztrsv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  100 INFOT = 1
      ztbsv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztbsv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztbsv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztbsv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      ztbsv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ztbsv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      ztbsv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  110 INFOT = 1
      ztpsv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztpsv('U', '/', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztpsv('U', 'N', '/', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztpsv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ztpsv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  120 INFOT = 1
      zgerc(-1, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgerc(0, -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgerc(0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zgerc(0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgerc(2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  130 INFOT = 1
      zgeru(-1, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgeru(0, -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgeru(0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zgeru(0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgeru(2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  140 INFOT = 1
      zher('/', 0, RALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zher('U', -1, RALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zher('U', 0, RALPHA, X, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zher('U', 2, RALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  150 INFOT = 1
      zhpr('/', 0, RALPHA, X, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhpr('U', -1, RALPHA, X, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zhpr('U', 0, RALPHA, X, 0, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  160 INFOT = 1
      zher2('/', 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zher2('U', -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zher2('U', 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zher2('U', 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      zher2('U', 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 180
  170 INFOT = 1
      zhpr2('/', 0, ALPHA, X, 1, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhpr2('U', -1, ALPHA, X, 1, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      zhpr2('U', 0, ALPHA, X, 0, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      zhpr2('U', 0, ALPHA, X, 1, Y, 0, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );

  180 IF( OK )THEN
         WRITE( NOUT, FMT = 9999 )SRNAMT
      } else {
         WRITE( NOUT, FMT = 9998 )SRNAMT
      }
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE TESTS OF ERROR-EXITS' )
 9998 FORMAT( ' ******* ', A6, ' FAILED THE TESTS OF ERROR-EXITS *****', '**' )

      // End of ZCHKE

      }
      SUBROUTINE ZMAKE( TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL )

*  Generates values for an M by N matrix A within the bandwidth
*  defined by KL and KU.
*  Stores the values in the array AA in the data structure required
*  by the routine, with unwanted elements set to rogue value.

*  TYPE is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D0, 0.0D0 ), ONE = ( 1.0D0, 0.0D0 ) ;
      COMPLEX*16         ROGUE
      const              ROGUE = ( -1.0D10, 1.0D10 ) ;
      double             RZERO;
      const              RZERO = 0.0D0 ;
      double             RROGUE;
      const              RROGUE = -1.0D10 ;
      // .. Scalar Arguments ..
      COMPLEX*16         TRANSL
      int                KL, KU, LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, * ), AA( * )
      // .. Local Scalars ..
      int                I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      COMPLEX*16         ZBEG
      // EXTERNAL ZBEG
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DCONJG, MAX, MIN
      // .. Executable Statements ..
      GEN = TYPE( 1: 1 ).EQ.'G'
      SYM = TYPE( 1: 1 ).EQ.'H'
      TRI = TYPE( 1: 1 ).EQ.'T'
      UPPER = ( SYM.OR.TRI ).AND.UPLO.EQ.'U'
      LOWER = ( SYM.OR.TRI ).AND.UPLO.EQ.'L'
      UNIT = TRI.AND.DIAG.EQ.'U'

      // Generate data in array A.

      DO 20 J = 1, N
         DO 10 I = 1, M
            if ( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) ) THEN                IF( ( I.LE.J.AND.J - I.LE.KU ).OR. ( I.GE.J.AND.I - J.LE.KL ) ) {
                  A( I, J ) = ZBEG( RESET ) + TRANSL
               } else {
                  A( I, J ) = ZERO
               }
               if ( I.NE.J ) {
                  if ( SYM ) {
                     A( J, I ) = DCONJG( A( I, J ) )
                  } else if ( TRI ) {
                     A( J, I ) = ZERO
                  }
               }
            }
   10    CONTINUE
         IF( SYM ) A( J, J ) = DCMPLX( DBLE( A( J, J ) ), RZERO )          IF( TRI ) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE
   20 CONTINUE

      // Store elements in array AS in data structure required by routine.

      if ( TYPE.EQ.'GE' ) {
         DO 50 J = 1, N
            DO 30 I = 1, M
               AA( I + ( J - 1 )*LDA ) = A( I, J )
   30       CONTINUE
            DO 40 I = M + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
   40       CONTINUE
   50    CONTINUE
      } else if ( TYPE.EQ.'GB' ) {
         DO 90 J = 1, N
            DO 60 I1 = 1, KU + 1 - J
               AA( I1 + ( J - 1 )*LDA ) = ROGUE
   60       CONTINUE
            DO 70 I2 = I1, MIN( KL + KU + 1, KU + 1 + M - J )
               AA( I2 + ( J - 1 )*LDA ) = A( I2 + J - KU - 1, J )
   70       CONTINUE
            DO 80 I3 = I2, LDA
               AA( I3 + ( J - 1 )*LDA ) = ROGUE
   80       CONTINUE
   90    CONTINUE
      } else if ( TYPE.EQ.'HE'.OR.TYPE.EQ.'TR' ) {
         DO 130 J = 1, N
            if ( UPPER ) {
               IBEG = 1
               if ( UNIT ) {
                  IEND = J - 1
               } else {
                  IEND = J
               }
            } else {
               if ( UNIT ) {
                  IBEG = J + 1
               } else {
                  IBEG = J
               }
               IEND = N
            }
            DO 100 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
  100       CONTINUE
            DO 110 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I, J )
  110       CONTINUE
            DO 120 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
  120       CONTINUE
            if ( SYM ) {
               JJ = J + ( J - 1 )*LDA
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE )
            }
  130    CONTINUE
      } else if ( TYPE.EQ.'HB'.OR.TYPE.EQ.'TB' ) {
         DO 170 J = 1, N
            if ( UPPER ) {
               KK = KL + 1
               IBEG = MAX( 1, KL + 2 - J )
               if ( UNIT ) {
                  IEND = KL
               } else {
                  IEND = KL + 1
               }
            } else {
               KK = 1
               if ( UNIT ) {
                  IBEG = 2
               } else {
                  IBEG = 1
               }
               IEND = MIN( KL + 1, 1 + M - J )
            }
            DO 140 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
  140       CONTINUE
            DO 150 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I + J - KK, J )
  150       CONTINUE
            DO 160 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
  160       CONTINUE
            if ( SYM ) {
               JJ = KK + ( J - 1 )*LDA
               AA( JJ ) = DCMPLX( DBLE( AA( JJ ) ), RROGUE )
            }
  170    CONTINUE
      } else if ( TYPE.EQ.'HP'.OR.TYPE.EQ.'TP' ) {
         IOFF = 0
         DO 190 J = 1, N
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            DO 180 I = IBEG, IEND
               IOFF = IOFF + 1
               AA( IOFF ) = A( I, J )
               if ( I.EQ.J ) {
                  IF( UNIT ) AA( IOFF ) = ROGUE                   IF( SYM ) AA( IOFF ) = DCMPLX( DBLE( AA( IOFF ) ), RROGUE )
               }
  180       CONTINUE
  190    CONTINUE
      }
      RETURN

      // End of ZMAKE

      }
      SUBROUTINE ZMVCH( TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, MV )

*  Checks the results of the computational tests.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D0, 0.0D0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0D0, RONE = 1.0D0 ;
      // .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      double             EPS, ERR;
      int                INCX, INCY, M, N, NMAX, NOUT;
      bool               FATAL, MV;
      String             TRANS;
      // .. Array Arguments ..
      COMPLEX*16         A( NMAX, * ), X( * ), Y( * ), YT( * ), YY( * )
      double             G( * );
      // .. Local Scalars ..
      COMPLEX*16         C
      double             ERRI;
      int                I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
      bool               CTRAN, TRAN;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, SQRT
      // .. Statement Functions ..
      double             ABS1;
      // .. Statement Function definitions ..
      ABS1( C ) = ABS( DBLE( C ) ) + ABS( DIMAG( C ) )
      // .. Executable Statements ..
      TRAN = TRANS.EQ.'T'
      CTRAN = TRANS.EQ.'C'
      if ( TRAN.OR.CTRAN ) {
         ML = N
         NL = M
      } else {
         ML = M
         NL = N
      }
      if ( INCX.LT.0 ) {
         KX = NL
         INCXL = -1
      } else {
         KX = 1
         INCXL = 1
      }
      if ( INCY.LT.0 ) {
         KY = ML
         INCYL = -1
      } else {
         KY = 1
         INCYL = 1
      }

      // Compute expected result in YT using data in A, X and Y.
      // Compute gauges in G.

      IY = KY
      DO 40 I = 1, ML
         YT( IY ) = ZERO
         G( IY ) = RZERO
         JX = KX
         if ( TRAN ) {
            DO 10 J = 1, NL
               YT( IY ) = YT( IY ) + A( J, I )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   10       CONTINUE
         } else if ( CTRAN ) {
            DO 20 J = 1, NL
               YT( IY ) = YT( IY ) + DCONJG( A( J, I ) )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   20       CONTINUE
         } else {
            DO 30 J = 1, NL
               YT( IY ) = YT( IY ) + A( I, J )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( I, J ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   30       CONTINUE
         }
         YT( IY ) = ALPHA*YT( IY ) + BETA*Y( IY )
         G( IY ) = ABS1( ALPHA )*G( IY ) + ABS1( BETA )*ABS1( Y( IY ) )
         IY = IY + INCYL
   40 CONTINUE

      // Compute the error ratio for this result.

      ERR = ZERO
      DO 50 I = 1, ML
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*ABS( INCY ) ) )/EPS
         IF( G( I ).NE.RZERO ) ERRI = ERRI/G( I )
         ERR = MAX( ERR, ERRI )
         IF( ERR*SQRT( EPS ).GE.RONE ) GO TO 60
   50 CONTINUE
      // If the loop completes, all results are at least half accurate.
      GO TO 80

      // Report fatal error.

   60 FATAL = .TRUE.
      WRITE( NOUT, FMT = 9999 )
      DO 70 I = 1, ML
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, YT( I ), YY( 1 + ( I - 1 )*ABS( INCY ) )
         } else {
            WRITE( NOUT, FMT = 9998 )I, YY( 1 + ( I - 1 )*ABS( INCY ) ), YT( I )
         }
   70 CONTINUE

   80 CONTINUE
      RETURN

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'                       EXPECTED RE', 'SULT                    COMPUTED RESULT' )
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) )

      // End of ZMVCH

      }
      bool    FUNCTION LZE( RI, RJ, LR );

*  Tests if two arrays are identical.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      COMPLEX*16         RI( * ), RJ( * )
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      DO 10 I = 1, LR
         IF( RI( I ).NE.RJ( I ) ) GO TO 20
   10 CONTINUE
      LZE = .TRUE.
      GO TO 30
   20 CONTINUE
      LZE = .FALSE.
   30 RETURN

      // End of LZE

      }
      bool    FUNCTION LZERES( TYPE, UPLO, M, N, AA, AS, LDA );

*  Tests if selected elements in two arrays are equal.

*  TYPE is 'GE', 'HE' or 'HP'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      COMPLEX*16         AA( LDA, * ), AS( LDA, * )
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO.EQ.'U'
      if ( TYPE.EQ.'GE' ) {
         DO 20 J = 1, N
            DO 10 I = M + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   10       CONTINUE
   20    CONTINUE
      } else if ( TYPE.EQ.'HE' ) {
         DO 50 J = 1, N
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            DO 30 I = 1, IBEG - 1
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   30       CONTINUE
            DO 40 I = IEND + 1, LDA
               IF( AA( I, J ).NE.AS( I, J ) ) GO TO 70
   40       CONTINUE
   50    CONTINUE
      }

      LZERES = .TRUE.
      GO TO 80
   70 CONTINUE
      LZERES = .FALSE.
   80 RETURN

      // End of LZERES

      }
      COMPLEX*16 FUNCTION ZBEG( RESET )

*  Generates complex numbers as pairs of random numbers uniformly
*  distributed between -0.5 and 0.5.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, J, MI, MJ;
      // .. Save statement ..
      SAVE               I, IC, J, MI, MJ
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891
         MJ = 457
         I = 7
         J = 7
         IC = 0
         RESET = .FALSE.
      }

      // The sequence of values of I or J is bounded between 1 and 999.
      // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I or J = 4 or 8, the period will be 25.
      // If initial I or J = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I or J
      // in 6.

      IC = IC + 1
   10 I = I*MI
      J = J*MJ
      I = I - 1000*( I/1000 )
      J = J - 1000*( J/1000 )
      if ( IC.GE.5 ) {
         IC = 0
         GO TO 10
      }
      ZBEG = DCMPLX( ( I - 500 )/1001.0D0, ( J - 500 )/1001.0D0 )
      RETURN

      // End of ZBEG

      }
      double           FUNCTION DDIFF( X, Y );

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.

      // .. Scalar Arguments ..
      double             X, Y;
      // .. Executable Statements ..
      DDIFF = X - Y
      RETURN

      // End of DDIFF

      }
      SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK )

*  Tests whether XERBLA has detected an error when it should.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Executable Statements ..
      if ( .NOT.LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT, SRNAMT
         OK = .FALSE.
      }
      LERR = .FALSE.
      RETURN

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ', I2, ' NOT D', 'ETECTED BY ', A6, ' *****' )

      // End of CHKXER

      }
      SUBROUTINE ZREGR1( TRANS, M, N, LY, KL, KU, ALPHA, A, LDA, X, INCX, BETA, Y, INCY, YS )

*  Input initialization for regression test.

      // .. Scalar Arguments ..
      String             TRANS;
      int                LY, M, N, KL, KU, LDA, INCX, INCY;
      COMPLEX*16         ALPHA, BETA
      // .. Array Arguments ..
      COMPLEX*16         A(LDA,*), X(*), Y(*), YS(*)
      // .. Local Scalars ..
      int                I;
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // .. Executable Statements ..
      TRANS = 'T'
      M = 0
      N = 5
      KL = 0
      KU = 0
      ALPHA = DCMPLX( 1.0D0 )
      LDA = MAX( 1, M )
      INCX = 1
      BETA = DCMPLX( -0.7D0, -0.8D0 )
      INCY = 1
      LY = ABS( INCY )*N
      DO 10 I = 1, LY
         Y( I ) = DCMPLX( 42.0D0, DBLE( I ) )
         YS( I ) = Y( I )
  10  CONTINUE
      RETURN
      }
      SUBROUTINE XERBLA( SRNAME, INFO )

*  This is a special version of XERBLA to be used only as part of
*  the test program for testing error exits from the Level 2 BLAS
*  routines.

*  XERBLA  is an error handler for the Level 2 BLAS routines.

*  It is called by the Level 2 BLAS routines if an input parameter is
*  invalid.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                INFO;
      String             SRNAME;
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Common blocks ..
      COMMON             /INFOC/INFOT, NOUT, OK, LERR
      COMMON             /SRNAMC/SRNAMT
      // .. Executable Statements ..
      LERR = .TRUE.
      if ( INFO.NE.INFOT ) {
         if ( INFOT.NE.0 ) {
            WRITE( NOUT, FMT = 9999 )INFO, INFOT
         } else {
            WRITE( NOUT, FMT = 9997 )INFO
         }
         OK = .FALSE.
      }
      if ( SRNAME.NE.SRNAMT ) {
         WRITE( NOUT, FMT = 9998 )SRNAME, SRNAMT
         OK = .FALSE.
      }
      RETURN

 9999 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' INSTEAD', ' OF ', I2, ' *******' )
 9998 FORMAT( ' ******* XERBLA WAS CALLED WITH SRNAME = ', A6, ' INSTE', 'AD OF ', A6, ' *******' )
 9997 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' *******' )

      // End of XERBLA

      }
