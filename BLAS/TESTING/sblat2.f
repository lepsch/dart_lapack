void main() {
*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

*  =====================================================================

      // .. Parameters ..
      int                NIN;
      const              NIN = 5 ;
      int                NSUBS;
      const              NSUBS = 16 ;
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      int                NMAX, INCMAX;
      const              NMAX = 65, INCMAX = 2 ;
      int                NINMAX, NIDMAX, NKBMAX, NALMAX, NBEMAX;
      const              NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7 ;
      // .. Local Scalars ..
      REAL               EPS, ERR, THRESH
      int                I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB, NOUT, NTRA;
      bool               FATAL, LTESTT, REWI, SAME, SFATAL, TRACE, TSTERR;
      String             TRANS;
      String             SNAMET;
      String             SNAPS, SUMMRY;
      // .. Local Arrays ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALMAX ), AS( NMAX*NMAX ), BET( NBEMAX ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( 2*NMAX )
      int                IDIM( NIDMAX ), INC( NINMAX ), KB( NKBMAX );
      bool               LTEST( NSUBS );
      String             SNAMES( NSUBS );
      // .. External Functions ..
      REAL               SDIFF
      bool               LSE;
      // EXTERNAL SDIFF, LSE
      // .. External Subroutines ..
      // EXTERNAL SCHK1, SCHK2, SCHK3, SCHK4, SCHK5, SCHK6, SCHKE, SMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      String             SRNAMT;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // COMMON /SRNAMC/SRNAMT
      // .. Data statements ..
      DATA               SNAMES/'SGEMV ', 'SGBMV ', 'SSYMV ', 'SSBMV ', 'SSPMV ', 'STRMV ', 'STBMV ', 'STPMV ', 'STRSV ', 'STBSV ', 'STPSV ', 'SGER  ', 'SSYR  ', 'SSPR  ', 'SSYR2 ', 'SSPR2 '/
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
      for (I = 1; I <= NIDIM; I++) { // 10
         if ( IDIM( I ).LT.0.OR.IDIM( I ).GT.NMAX ) {
            WRITE( NOUT, FMT = 9996 )NMAX
            GO TO 230
         }
      } // 10
      // Values of K
      READ( NIN, FMT = * )NKB
      if ( NKB.LT.1.OR.NKB.GT.NKBMAX ) {
         WRITE( NOUT, FMT = 9997 )'K', NKBMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( KB( I ), I = 1, NKB )
      for (I = 1; I <= NKB; I++) { // 20
         if ( KB( I ).LT.0 ) {
            WRITE( NOUT, FMT = 9995 )
            GO TO 230
         }
      } // 20
      // Values of INCX and INCY
      READ( NIN, FMT = * )NINC
      if ( NINC.LT.1.OR.NINC.GT.NINMAX ) {
         WRITE( NOUT, FMT = 9997 )'INCX AND INCY', NINMAX
         GO TO 230
      }
      READ( NIN, FMT = * )( INC( I ), I = 1, NINC )
      for (I = 1; I <= NINC; I++) { // 30
         if ( INC( I ) == 0.OR.ABS( INC( I ) ).GT.INCMAX ) {
            WRITE( NOUT, FMT = 9994 )INCMAX
            GO TO 230
         }
      } // 30
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

      for (I = 1; I <= NSUBS; I++) { // 40
         LTEST( I ) = false;
      } // 40
   50 READ( NIN, FMT = 9984, END = 80 )SNAMET, LTESTT
      for (I = 1; I <= NSUBS; I++) { // 60
         IF( SNAMET == SNAMES( I ) ) GO TO 70
      } // 60
      WRITE( NOUT, FMT = 9986 )SNAMET
      STOP
   70 LTEST( I ) = LTESTT
      GO TO 50

      } // 80
      CLOSE ( NIN )

      // Compute EPS (the machine precision).

      EPS = EPSILON(ZERO)
      WRITE( NOUT, FMT = 9998 )EPS

      // Check the reliability of SMVCH using exact data.

      N = MIN( 32, NMAX )
      for (J = 1; J <= N; J++) { // 120
         for (I = 1; I <= N; I++) { // 110
            A( I, J ) = MAX( I - J + 1, 0 )
         } // 110
         X( J ) = J
         Y( J ) = ZERO
      } // 120
      for (J = 1; J <= N; J++) { // 130
         YY( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
      } // 130
      // YY holds the exact result. On exit from SMVCH YT holds
      // the result computed by SMVCH.
      TRANS = 'N'
      smvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR != ZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }
      TRANS = 'T'
      smvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
      SAME = LSE( YY, YT, N )
      if ( .NOT.SAME.OR.ERR != ZERO ) {
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
         STOP
      }

      // Test each subroutine in turn.

      for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) { // 210
         WRITE( NOUT, FMT = * )
         if ( .NOT.LTEST( ISNUM ) ) {
            // Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM )
         } else {
            SRNAMT = SNAMES( ISNUM )
            // Test error exits.
            if ( TSTERR ) {
               schke(ISNUM, SNAMES( ISNUM ), NOUT );
               WRITE( NOUT, FMT = * )
            }
            // Test computations.
            INFOT = 0
            OK = true;
            FATAL = false;
            GO TO ( 140, 140, 150, 150, 150, 160, 160, 160, 160, 160, 160, 170, 180, 180, 190, 190 )ISNUM
            // Test SGEMV, 01, and SGBMV, 02.
  140       CALL SCHK1( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
            // Test SSYMV, 03, SSBMV, 04, and SSPMV, 05.
  150       CALL SCHK2( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
            // Test STRMV, 06, STBMV, 07, STPMV, 08,
            // STRSV, 09, STBSV, 10, and STPSV, 11.
  160       CALL SCHK3( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test SGER, 12.
  170       CALL SCHK4( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test SSYR, 13, and SSPR, 14.
  180       CALL SCHK5( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )
            GO TO 200
            // Test SSYR2, 15, and SSPR2, 16.
  190       CALL SCHK6( SNAMES( ISNUM ), EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

  200       IF( FATAL.AND.SFATAL ) GO TO 220
         }
      } // 210
      WRITE( NOUT, FMT = 9982 )
      GO TO 240

      } // 220
      WRITE( NOUT, FMT = 9981 )
      GO TO 240

      } // 230
      WRITE( NOUT, FMT = 9987 )

      } // 240
      if (TRACE) CLOSE ( NTRA );
      CLOSE ( NOUT )
      STOP

 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES', 'S THAN', F8.2 )
 9998 FORMAT( ' RELATIVE MACHINE PRECISION IS TAKEN TO BE', 1P, E9.1 )
 9997 FORMAT( ' NUMBER OF VALUES OF ', A, ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( ' VALUE OF N IS LESS THAN 0 OR GREATER THAN ', I2 )
 9995 FORMAT( ' VALUE OF K IS LESS THAN 0' )
 9994 FORMAT( ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ', I2 )
 9993 FORMAT( ' TESTS OF THE REAL             LEVEL 2 BLAS', //' THE F', 'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9992 FORMAT( '   FOR N              ', 9I6 )
 9991 FORMAT( '   FOR K              ', 7I6 )
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 )
 9989 FORMAT( '   FOR ALPHA          ', 7F6.1 )
 9988 FORMAT( '   FOR BETA           ', 7F6.1 )
 9987 FORMAT( ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM', /' ******* TESTS ABANDONED *******' )
 9986 FORMAT( ' SUBPROGRAM NAME ', A6, ' NOT RECOGNIZED', /' ******* T', 'ESTS ABANDONED *******' )
 9985 FORMAT( ' ERROR IN SMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU', 'ATED WRONGLY.', /' SMVCH WAS CALLED WITH TRANS = ', A1, ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '.', / ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.' , /' ******* TESTS ABANDONED *******' )
 9984 FORMAT( A6, L2 )
 9983 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9982 FORMAT( /' END OF TESTS' )
 9981 FORMAT( /' ******* FATAL ERROR - TESTS ABANDONED *******' )
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )

      // End of SBLAT2

      }
      SUBROUTINE SCHK1( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )

*  Tests SGEMV and SGBMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF
      const              ZERO = 0.0, HALF = 0.5 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      REAL               ALPHA, ALS, BETA, BLS, ERR, ERRMAX, TRANSL
      int                I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA, LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK, NL, NS;
      bool               BANDED, FULL, NULL, RESET, SAME, TRAN;
      String             TRANS, TRANSS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SGBMV, SGEMV, SMAKE, SMVCH, SREGR1
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'NTC'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'E'
      BANDED = SNAME( 3: 3 ) == 'B'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 11
      } else if ( BANDED ) {
         NARGS = 13
      }

      NC = 0
      RESET = true;
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN )
         ND = N/2 + 1

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = MAX( N - ND, 0 )             IF( IM == 2 ) M = MIN( N + ND, NMAX );

            if ( BANDED ) {
               NK = NKB
            } else {
               NK = 1
            }
            for (IKU = 1; IKU <= NK; IKU++) { // 100
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
               if (LDA.LT.NMAX) LDA = LDA + 1;
               // Skip tests if not enough room.
               if (LDA.GT.NMAX) GO TO 100;
               LAA = LDA*N
               NULL = N.LE.0.OR.M.LE.0

               // Generate the matrix A.

               TRANSL = ZERO
               smake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL );

               for (IC = 1; IC <= 3; IC++) { // 90
                  TRANS = ICH( IC: IC )
                  TRAN = TRANS == 'T'.OR.TRANS == 'C'

                  if ( TRAN ) {
                     ML = N
                     NL = M
                  } else {
                     ML = M
                     NL = N
                  }

                  for (IX = 1; IX <= NINC; IX++) { // 80
                     INCX = INC( IX )
                     LX = ABS( INCX )*NL

                     // Generate the vector X.

                     TRANSL = HALF
                     smake('GE', ' ', ' ', 1, NL, X, 1, XX, ABS( INCX ), 0, NL - 1, RESET, TRANSL );
                     if ( NL.GT.1 ) {
                        X( NL/2 ) = ZERO
                        XX( 1 + ABS( INCX )*( NL/2 - 1 ) ) = ZERO
                     }

                     for (IY = 1; IY <= NINC; IY++) { // 70
                        INCY = INC( IY )
                        LY = ABS( INCY )*ML

                        for (IA = 1; IA <= NALF; IA++) { // 60
                           ALPHA = ALF( IA )

                           for (IB = 1; IB <= NBET; IB++) { // 50
                              BETA = BET( IB )

                              // Generate the vector Y.

                              TRANSL = ZERO
                              smake('GE', ' ', ' ', 1, ML, Y, 1, YY, ABS( INCY ), 0, ML - 1, RESET, TRANSL );

                              NC = NC + 1

                              // Save every datum before calling the
                              // subroutine.

                              TRANSS = TRANS
                              MS = M
                              NS = N
                              KLS = KL
                              KUS = KU
                              ALS = ALPHA
                              for (I = 1; I <= LAA; I++) { // 10
                                 AS( I ) = AA( I )
                              } // 10
                              LDAS = LDA
                              for (I = 1; I <= LX; I++) { // 20
                                 XS( I ) = XX( I )
                              } // 20
                              INCXS = INCX
                              BLS = BETA
                              for (I = 1; I <= LY; I++) { // 30
                                 YS( I ) = YY( I )
                              } // 30
                              INCYS = INCY

                              // Call the subroutine.

                              if ( FULL ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 sgemv(TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              } else if ( BANDED ) {
                                 if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
                                 if (REWI) REWIND NTRA;
                                 sgbmv(TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                              }

                              // Check if error-exit was taken incorrectly.

                              if ( .NOT.OK ) {
                                 WRITE( NOUT, FMT = 9993 )
                                 FATAL = true;
                                 GO TO 130
                              }

                              // See what data changed inside subroutines.

                              ISAME( 1 ) = TRANS == TRANSS
                              ISAME( 2 ) = MS == M
                              ISAME( 3 ) = NS == N
                              if ( FULL ) {
                                 ISAME( 4 ) = ALS == ALPHA
                                 ISAME( 5 ) = LSE( AS, AA, LAA )
                                 ISAME( 6 ) = LDAS == LDA
                                 ISAME( 7 ) = LSE( XS, XX, LX )
                                 ISAME( 8 ) = INCXS == INCX
                                 ISAME( 9 ) = BLS == BETA
                                 if ( NULL ) {
                                    ISAME( 10 ) = LSE( YS, YY, LY )
                                 } else {
                                    ISAME( 10 ) = LSERES( 'GE', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 11 ) = INCYS == INCY
                              } else if ( BANDED ) {
                                 ISAME( 4 ) = KLS == KL
                                 ISAME( 5 ) = KUS == KU
                                 ISAME( 6 ) = ALS == ALPHA
                                 ISAME( 7 ) = LSE( AS, AA, LAA )
                                 ISAME( 8 ) = LDAS == LDA
                                 ISAME( 9 ) = LSE( XS, XX, LX )
                                 ISAME( 10 ) = INCXS == INCX
                                 ISAME( 11 ) = BLS == BETA
                                 if ( NULL ) {
                                    ISAME( 12 ) = LSE( YS, YY, LY )
                                 } else {
                                    ISAME( 12 ) = LSERES( 'GE', ' ', 1, ML, YS, YY, ABS( INCY ) )
                                 }
                                 ISAME( 13 ) = INCYS == INCY
                              }

                              // If data was incorrectly changed, report
                              // and return.

                              SAME = true;
                              for (I = 1; I <= NARGS; I++) { // 40
                                 SAME = SAME.AND.ISAME( I )
                                 IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                              } // 40
                              if ( .NOT.SAME ) {
                                 FATAL = true;
                                 GO TO 130
                              }

                              if ( .NOT.NULL ) {

                                 // Check the result.

                                 smvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                                 ERRMAX = MAX( ERRMAX, ERR )
                                 // If got really bad answer, report and
                                 // return.
                                 if (FATAL) GO TO 130;
                              } else {
                                 // Avoid repeating tests with M.le.0 or
                                 // N.le.0.
                                 GO TO 110
                              }

                           } // 50

                        } // 60

                     } // 70

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Regression test to verify preservation of y when m zero, n nonzero.

      sregr1(TRANS, M, N, LY, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY, YS );
      if ( FULL ) {
         if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY;
         if (REWI) REWIND NTRA;
         sgemv(TRANS, M, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
      } else if ( BANDED ) {
         if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY;
         if (REWI) REWIND NTRA;
         sgbmv(TRANS, M, N, KL, KU, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
      }
      NC = NC + 1
      if ( .NOT.LSE( YS, YY, LY ) ) {
         WRITE( NOUT, FMT = 9998 )NARGS - 1
         FATAL = true;
         GO TO 130
      }

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 140

      } // 130
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY
      }

      } // 140
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 4( I3, ',' ), F4.1, ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2, ') .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), F4.1, ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2, ')         .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK1

      }
      SUBROUTINE SCHK2( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G )

*  Tests SSYMV, SSBMV and SSPMV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF
      const              ZERO = 0.0, HALF = 0.5 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NALF, NBET, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX )
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      REAL               ALPHA, ALS, BETA, BLS, ERR, ERRMAX, TRANSL
      int                I, IA, IB, IC, IK, IN, INCX, INCXS, INCY, INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMVCH, SSBMV, SSPMV, SSYMV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'Y'
      BANDED = SNAME( 3: 3 ) == 'B'
      PACKED = SNAME( 3: 3 ) == 'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 10
      } else if ( BANDED ) {
         NARGS = 11
      } else if ( PACKED ) {
         NARGS = 9
      }

      NC = 0
      RESET = true;
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         for (IK = 1; IK <= NK; IK++) { // 100
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
            if (LDA.LT.NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA.GT.NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            for (IC = 1; IC <= 2; IC++) { // 90
               UPLO = ICH( IC: IC )

               // Generate the matrix A.

               TRANSL = ZERO
               smake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

               for (IX = 1; IX <= NINC; IX++) { // 80
                  INCX = INC( IX )
                  LX = ABS( INCX )*N

                  // Generate the vector X.

                  TRANSL = HALF
                  smake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     X( N/2 ) = ZERO
                     XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IY = 1; IY <= NINC; IY++) { // 70
                     INCY = INC( IY )
                     LY = ABS( INCY )*N

                     for (IA = 1; IA <= NALF; IA++) { // 60
                        ALPHA = ALF( IA )

                        for (IB = 1; IB <= NBET; IB++) { // 50
                           BETA = BET( IB )

                           // Generate the vector Y.

                           TRANSL = ZERO
                           smake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );

                           NC = NC + 1

                           // Save every datum before calling the
                           // subroutine.

                           UPLOS = UPLO
                           NS = N
                           KS = K
                           ALS = ALPHA
                           for (I = 1; I <= LAA; I++) { // 10
                              AS( I ) = AA( I )
                           } // 10
                           LDAS = LDA
                           for (I = 1; I <= LX; I++) { // 20
                              XS( I ) = XX( I )
                           } // 20
                           INCXS = INCX
                           BLS = BETA
                           for (I = 1; I <= LY; I++) { // 30
                              YS( I ) = YY( I )
                           } // 30
                           INCYS = INCY

                           // Call the subroutine.

                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              ssymv(UPLO, N, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              ssbmv(UPLO, N, K, ALPHA, AA, LDA, XX, INCX, BETA, YY, INCY );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY;
                              if (REWI) REWIND NTRA;
                              sspmv(UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY );
                           }

                           // Check if error-exit was taken incorrectly.

                           if ( .NOT.OK ) {
                              WRITE( NOUT, FMT = 9992 )
                              FATAL = true;
                              GO TO 120
                           }

                           // See what data changed inside subroutines.

                           ISAME( 1 ) = UPLO == UPLOS
                           ISAME( 2 ) = NS == N
                           if ( FULL ) {
                              ISAME( 3 ) = ALS == ALPHA
                              ISAME( 4 ) = LSE( AS, AA, LAA )
                              ISAME( 5 ) = LDAS == LDA
                              ISAME( 6 ) = LSE( XS, XX, LX )
                              ISAME( 7 ) = INCXS == INCX
                              ISAME( 8 ) = BLS == BETA
                              if ( NULL ) {
                                 ISAME( 9 ) = LSE( YS, YY, LY )
                              } else {
                                 ISAME( 9 ) = LSERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 10 ) = INCYS == INCY
                           } else if ( BANDED ) {
                              ISAME( 3 ) = KS == K
                              ISAME( 4 ) = ALS == ALPHA
                              ISAME( 5 ) = LSE( AS, AA, LAA )
                              ISAME( 6 ) = LDAS == LDA
                              ISAME( 7 ) = LSE( XS, XX, LX )
                              ISAME( 8 ) = INCXS == INCX
                              ISAME( 9 ) = BLS == BETA
                              if ( NULL ) {
                                 ISAME( 10 ) = LSE( YS, YY, LY )
                              } else {
                                 ISAME( 10 ) = LSERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 11 ) = INCYS == INCY
                           } else if ( PACKED ) {
                              ISAME( 3 ) = ALS == ALPHA
                              ISAME( 4 ) = LSE( AS, AA, LAA )
                              ISAME( 5 ) = LSE( XS, XX, LX )
                              ISAME( 6 ) = INCXS == INCX
                              ISAME( 7 ) = BLS == BETA
                              if ( NULL ) {
                                 ISAME( 8 ) = LSE( YS, YY, LY )
                              } else {
                                 ISAME( 8 ) = LSERES( 'GE', ' ', 1, N, YS, YY, ABS( INCY ) )
                              }
                              ISAME( 9 ) = INCYS == INCY
                           }

                           // If data was incorrectly changed, report and
                           // return.

                           SAME = true;
                           for (I = 1; I <= NARGS; I++) { // 40
                              SAME = SAME.AND.ISAME( I )
                              IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                           } // 40
                           if ( .NOT.SAME ) {
                              FATAL = true;
                              GO TO 120
                           }

                           if ( .NOT.NULL ) {

                              // Check the result.

                              smvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, true );
                              ERRMAX = MAX( ERRMAX, ERR )
                              // If got really bad answer, report and
                              // return.
                              if (FATAL) GO TO 120;
                           } else {
                              // Avoid repeating tests with N.le.0
                              GO TO 110
                           }

                        } // 50

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY
      }

      } // 130
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', AP', ', X,', I2, ',', F4.1, ', Y,', I2, ')                .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), F4.1, ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2, ')         .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2, ')             .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK2

      }
      SUBROUTINE SCHK3( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, XT, G, Z )

*  Tests STRMV, STBMV, STPMV, STRSV, STBSV and STPSV.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NIDIM, NINC, NKB, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), AS( NMAX*NMAX ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XT( NMAX ), XX( NMAX*INCMAX ), Z( NMAX )
      int                IDIM( NIDIM ), INC( NINC ), KB( NKB );
      // .. Local Scalars ..
      REAL               ERR, ERRMAX, TRANSL
      int                I, ICD, ICT, ICU, IK, IN, INCX, INCXS, IX, K, KS, LAA, LDA, LDAS, LX, N, NARGS, NC, NK, NS;
      bool               BANDED, FULL, NULL, PACKED, RESET, SAME;
      String             DIAG, DIAGS, TRANS, TRANSS, UPLO, UPLOS;
      String             ICHD, ICHU;
      String             ICHT;
      // .. Local Arrays ..
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMVCH, STBMV, STBSV, STPMV, STPSV, STRMV, STRSV
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'R'
      BANDED = SNAME( 3: 3 ) == 'B'
      PACKED = SNAME( 3: 3 ) == 'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 8
      } else if ( BANDED ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 7
      }

      NC = 0
      RESET = true;
      ERRMAX = ZERO
      // Set up zero vector for SMVCH.
      for (I = 1; I <= NMAX; I++) { // 10
         Z( I ) = ZERO
      } // 10

      for (IN = 1; IN <= NIDIM; IN++) { // 110
         N = IDIM( IN )

         if ( BANDED ) {
            NK = NKB
         } else {
            NK = 1
         }
         for (IK = 1; IK <= NK; IK++) { // 100
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
            if (LDA.LT.NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA.GT.NMAX) GO TO 100;
            if ( PACKED ) {
               LAA = ( N*( N + 1 ) )/2
            } else {
               LAA = LDA*N
            }
            NULL = N.LE.0

            for (ICU = 1; ICU <= 2; ICU++) { // 90
               UPLO = ICHU( ICU: ICU )

               for (ICT = 1; ICT <= 3; ICT++) { // 80
                  TRANS = ICHT( ICT: ICT )

                  for (ICD = 1; ICD <= 2; ICD++) { // 70
                     DIAG = ICHD( ICD: ICD )

                     // Generate the matrix A.

                     TRANSL = ZERO
                     smake(SNAME( 2: 3 ), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K, K, RESET, TRANSL );

                     for (IX = 1; IX <= NINC; IX++) { // 60
                        INCX = INC( IX )
                        LX = ABS( INCX )*N

                        // Generate the vector X.

                        TRANSL = HALF
                        smake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
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
                        for (I = 1; I <= LAA; I++) { // 20
                           AS( I ) = AA( I )
                        } // 20
                        LDAS = LDA
                        for (I = 1; I <= LX; I++) { // 30
                           XS( I ) = XX( I )
                        } // 30
                        INCXS = INCX

                        // Call the subroutine.

                        if ( SNAME( 4: 5 ) == 'MV' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              strmv(UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              stbmv(UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              stpmv(UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        } else if ( SNAME( 4: 5 ) == 'SV' ) {
                           if ( FULL ) {
                              if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              strsv(UPLO, TRANS, DIAG, N, AA, LDA, XX, INCX );
                           } else if ( BANDED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX;
                              if (REWI) REWIND NTRA;
                              stbsv(UPLO, TRANS, DIAG, N, K, AA, LDA, XX, INCX );
                           } else if ( PACKED ) {
                              if (TRACE) WRITE( NTRA, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX;
                              if (REWI) REWIND NTRA;
                              stpsv(UPLO, TRANS, DIAG, N, AA, XX, INCX );
                           }
                        }

                        // Check if error-exit was taken incorrectly.

                        if ( .NOT.OK ) {
                           WRITE( NOUT, FMT = 9992 )
                           FATAL = true;
                           GO TO 120
                        }

                        // See what data changed inside subroutines.

                        ISAME( 1 ) = UPLO == UPLOS
                        ISAME( 2 ) = TRANS == TRANSS
                        ISAME( 3 ) = DIAG == DIAGS
                        ISAME( 4 ) = NS == N
                        if ( FULL ) {
                           ISAME( 5 ) = LSE( AS, AA, LAA )
                           ISAME( 6 ) = LDAS == LDA
                           if ( NULL ) {
                              ISAME( 7 ) = LSE( XS, XX, LX )
                           } else {
                              ISAME( 7 ) = LSERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 8 ) = INCXS == INCX
                        } else if ( BANDED ) {
                           ISAME( 5 ) = KS == K
                           ISAME( 6 ) = LSE( AS, AA, LAA )
                           ISAME( 7 ) = LDAS == LDA
                           if ( NULL ) {
                              ISAME( 8 ) = LSE( XS, XX, LX )
                           } else {
                              ISAME( 8 ) = LSERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 9 ) = INCXS == INCX
                        } else if ( PACKED ) {
                           ISAME( 5 ) = LSE( AS, AA, LAA )
                           if ( NULL ) {
                              ISAME( 6 ) = LSE( XS, XX, LX )
                           } else {
                              ISAME( 6 ) = LSERES( 'GE', ' ', 1, N, XS, XX, ABS( INCX ) )
                           }
                           ISAME( 7 ) = INCXS == INCX
                        }

                        // If data was incorrectly changed, report and
                        // return.

                        SAME = true;
                        for (I = 1; I <= NARGS; I++) { // 40
                           SAME = SAME.AND.ISAME( I )
                           IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                        } // 40
                        if ( .NOT.SAME ) {
                           FATAL = true;
                           GO TO 120
                        }

                        if ( .NOT.NULL ) {
                           if ( SNAME( 4: 5 ) == 'MV' ) {

                              // Check the result.

                              smvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, true );
                           } else if ( SNAME( 4: 5 ) == 'SV' ) {

                              // Compute approximation to original vector.

                              for (I = 1; I <= N; I++) { // 50
                                 Z( I ) = XX( 1 + ( I - 1 )* ABS( INCX ) )                                  XX( 1 + ( I - 1 )*ABS( INCX ) ) = X( I )
                              } // 50
                              smvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT, G, XX, EPS, ERR, FATAL, NOUT, false );
                           }
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 120;
                        } else {
                           // Avoid repeating tests with N.le.0.
                           GO TO 110
                        }

                     } // 60

                  } // 70

               } // 80

            } // 90

         } // 100

      } // 110

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX
      } else if ( BANDED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9995 )NC, SNAME, UPLO, TRANS, DIAG, N, INCX
      }

      } // 130
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', AP, ', 'X,', I2, ')                        .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), 2( I3, ',' ), ' A,', I3, ', X,', I2, ')                 .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', A,', I3, ', X,', I2, ')                     .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK3

      }
      SUBROUTINE SCHK4( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests SGER.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      REAL               ALPHA, ALS, ERR, ERRMAX, TRANSL
      int                I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS, NC, ND, NS;
      bool               NULL, RESET, SAME;
      // .. Local Arrays ..
      REAL               W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SGER, SMAKE, SMVCH
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      // Define the number of arguments.
      NARGS = 9

      NC = 0
      RESET = true;
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 120
         N = IDIM( IN )
         ND = N/2 + 1

         for (IM = 1; IM <= 2; IM++) { // 110
            if (IM == 1) M = MAX( N - ND, 0 )             IF( IM == 2 ) M = MIN( N + ND, NMAX );

            // Set LDA to 1 more than minimum value if room.
            LDA = M
            if (LDA.LT.NMAX) LDA = LDA + 1;
            // Skip tests if not enough room.
            if (LDA.GT.NMAX) GO TO 110;
            LAA = LDA*N
            NULL = N.LE.0.OR.M.LE.0

            for (IX = 1; IX <= NINC; IX++) { // 100
               INCX = INC( IX )
               LX = ABS( INCX )*M

               // Generate the vector X.

               TRANSL = HALF
               smake('GE', ' ', ' ', 1, M, X, 1, XX, ABS( INCX ), 0, M - 1, RESET, TRANSL );
               if ( M.GT.1 ) {
                  X( M/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( M/2 - 1 ) ) = ZERO
               }

               for (IY = 1; IY <= NINC; IY++) { // 90
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  smake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 80
                     ALPHA = ALF( IA )

                     // Generate the matrix A.

                     TRANSL = ZERO
                     smake(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA, LDA, M - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     MS = M
                     NS = N
                     ALS = ALPHA
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I )
                     } // 10
                     LDAS = LDA
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I )
                     } // 20
                     INCXS = INCX
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I )
                     } // 30
                     INCYS = INCY

                     // Call the subroutine.

                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA;
                     if (REWI) REWIND NTRA;
                     sger(M, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9993 )
                        FATAL = true;
                        GO TO 140
                     }

                     // See what data changed inside subroutine.

                     ISAME( 1 ) = MS == M
                     ISAME( 2 ) = NS == N
                     ISAME( 3 ) = ALS == ALPHA
                     ISAME( 4 ) = LSE( XS, XX, LX )
                     ISAME( 5 ) = INCXS == INCX
                     ISAME( 6 ) = LSE( YS, YY, LY )
                     ISAME( 7 ) = INCYS == INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LSE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LSERES( 'GE', ' ', M, N, AS, AA, LDA )
                     }
                     ISAME( 9 ) = LDAS == LDA

                     // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                     } // 40
                     if ( .NOT.SAME ) {
                        FATAL = true;
                        GO TO 140
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           for (I = 1; I <= M; I++) { // 50
                              Z( I ) = X( I )
                           } // 50
                        } else {
                           for (I = 1; I <= M; I++) { // 60
                              Z( I ) = X( M - I + 1 )
                           } // 60
                        }
                        for (J = 1; J <= N; J++) { // 70
                           if ( INCY.GT.0 ) {
                              W( 1 ) = Y( J )
                           } else {
                              W( 1 ) = Y( N - J + 1 )
                           }
                           smvch('N', M, 1, ALPHA, Z, NMAX, W, 1, ONE, A( 1, J ), 1, YT, G, AA( 1 + ( J - 1 )*LDA ), EPS, ERR, FATAL, NOUT, true );
                           ERRMAX = MAX( ERRMAX, ERR )
                           // If got really bad answer, report and return.
                           if (FATAL) GO TO 130;
                        } // 70
                     } else {
                        // Avoid repeating tests with M.le.0 or N.le.0.
                        GO TO 110
                     }

                  } // 80

               } // 90

            } // 100

         } // 110

      } // 120

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 150

      } // 130
      WRITE( NOUT, FMT = 9995 )J

      } // 140
      WRITE( NOUT, FMT = 9996 )SNAME
      WRITE( NOUT, FMT = 9994 )NC, SNAME, M, N, ALPHA, INCX, INCY, LDA

      } // 150
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( I3, ',' ), F4.1, ', X,', I2, ', Y,', I2, ', A,', I3, ')                  .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK4

      }
      SUBROUTINE SCHK5( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests SSYR and SSPR.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX )
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      REAL               ALPHA, ALS, ERR, ERRMAX, TRANSL
      int                I, IA, IC, IN, INCX, INCXS, IX, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      REAL               W( 1 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMVCH, SSPR, SSYR
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'Y'
      PACKED = SNAME( 3: 3 ) == 'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 7
      } else if ( PACKED ) {
         NARGS = 6
      }

      NC = 0
      RESET = true;
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 100
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         if (LDA.LT.NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA.GT.NMAX) GO TO 100;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         for (IC = 1; IC <= 2; IC++) { // 90
            UPLO = ICH( IC: IC )
            UPPER = UPLO == 'U'

            for (IX = 1; IX <= NINC; IX++) { // 80
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               smake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               for (IA = 1; IA <= NALF; IA++) { // 70
                  ALPHA = ALF( IA )
                  NULL = N.LE.0.OR.ALPHA == ZERO

                  // Generate the matrix A.

                  TRANSL = ZERO
                  smake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                  NC = NC + 1

                  // Save every datum before calling the subroutine.

                  UPLOS = UPLO
                  NS = N
                  ALS = ALPHA
                  for (I = 1; I <= LAA; I++) { // 10
                     AS( I ) = AA( I )
                  } // 10
                  LDAS = LDA
                  for (I = 1; I <= LX; I++) { // 20
                     XS( I ) = XX( I )
                  } // 20
                  INCXS = INCX

                  // Call the subroutine.

                  if ( FULL ) {
                     if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, LDA;
                     if (REWI) REWIND NTRA;
                     ssyr(UPLO, N, ALPHA, XX, INCX, AA, LDA );
                  } else if ( PACKED ) {
                     if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX;
                     if (REWI) REWIND NTRA;
                     sspr(UPLO, N, ALPHA, XX, INCX, AA );
                  }

                  // Check if error-exit was taken incorrectly.

                  if ( .NOT.OK ) {
                     WRITE( NOUT, FMT = 9992 )
                     FATAL = true;
                     GO TO 120
                  }

                  // See what data changed inside subroutines.

                  ISAME( 1 ) = UPLO == UPLOS
                  ISAME( 2 ) = NS == N
                  ISAME( 3 ) = ALS == ALPHA
                  ISAME( 4 ) = LSE( XS, XX, LX )
                  ISAME( 5 ) = INCXS == INCX
                  if ( NULL ) {
                     ISAME( 6 ) = LSE( AS, AA, LAA )
                  } else {
                     ISAME( 6 ) = LSERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA )
                  }
                  if ( .NOT.PACKED ) {
                     ISAME( 7 ) = LDAS == LDA
                  }

                  // If data was incorrectly changed, report and return.

                  SAME = true;
                  for (I = 1; I <= NARGS; I++) { // 30
                     SAME = SAME.AND.ISAME( I )
                     IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                  } // 30
                  if ( .NOT.SAME ) {
                     FATAL = true;
                     GO TO 120
                  }

                  if ( .NOT.NULL ) {

                     // Check the result column by column.

                     if ( INCX.GT.0 ) {
                        for (I = 1; I <= N; I++) { // 40
                           Z( I ) = X( I )
                        } // 40
                     } else {
                        for (I = 1; I <= N; I++) { // 50
                           Z( I ) = X( N - I + 1 )
                        } // 50
                     }
                     JA = 1
                     for (J = 1; J <= N; J++) { // 60
                        W( 1 ) = Z( J )
                        if ( UPPER ) {
                           JJ = 1
                           LJ = J
                        } else {
                           JJ = J
                           LJ = N - J + 1
                        }
                        smvch('N', LJ, 1, ALPHA, Z( JJ ), LJ, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
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
                        if (FATAL) GO TO 110;
                     } // 60
                  } else {
                     // Avoid repeating tests if N.le.0.
                     if (N.LE.0) GO TO 100;
                  }

               } // 70

            } // 80

         } // 90

      } // 100

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 130

      } // 110
      WRITE( NOUT, FMT = 9995 )J

      } // 120
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX
      }

      } // 130
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', AP)                           .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', A,', I3, ')                        .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK5

      }
      SUBROUTINE SCHK6( SNAME, EPS, THRESH, NOUT, NTRA, TRACE, REWI, FATAL, NIDIM, IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS, YT, G, Z )

*  Tests SSYR2 and SSPR2.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
      // .. Scalar Arguments ..
      REAL               EPS, THRESH
      int                INCMAX, NALF, NIDIM, NINC, NMAX, NOUT, NTRA;
      bool               FATAL, REWI, TRACE;
      String             SNAME;
      // .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ), AS( NMAX*NMAX ), G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ), XX( NMAX*INCMAX ), Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ), YY( NMAX*INCMAX ), Z( NMAX, 2 )
      int                IDIM( NIDIM ), INC( NINC );
      // .. Local Scalars ..
      REAL               ALPHA, ALS, ERR, ERRMAX, TRANSL
      int                I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX, IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N, NARGS, NC, NS;
      bool               FULL, NULL, PACKED, RESET, SAME, UPPER;
      String             UPLO, UPLOS;
      String             ICH;
      // .. Local Arrays ..
      REAL               W( 2 )
      bool               ISAME( 13 );
      // .. External Functions ..
      bool               LSE, LSERES;
      // EXTERNAL LSE, LSERES
      // .. External Subroutines ..
      // EXTERNAL SMAKE, SMVCH, SSPR2, SSYR2
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // .. Scalars in Common ..
      int                INFOT, NOUTC;
      bool               LERR, OK;
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Data statements ..
      DATA               ICH/'UL'/
      // .. Executable Statements ..
      FULL = SNAME( 3: 3 ) == 'Y'
      PACKED = SNAME( 3: 3 ) == 'P'
      // Define the number of arguments.
      if ( FULL ) {
         NARGS = 9
      } else if ( PACKED ) {
         NARGS = 8
      }

      NC = 0
      RESET = true;
      ERRMAX = ZERO

      for (IN = 1; IN <= NIDIM; IN++) { // 140
         N = IDIM( IN )
         // Set LDA to 1 more than minimum value if room.
         LDA = N
         if (LDA.LT.NMAX) LDA = LDA + 1;
         // Skip tests if not enough room.
         if (LDA.GT.NMAX) GO TO 140;
         if ( PACKED ) {
            LAA = ( N*( N + 1 ) )/2
         } else {
            LAA = LDA*N
         }

         for (IC = 1; IC <= 2; IC++) { // 130
            UPLO = ICH( IC: IC )
            UPPER = UPLO == 'U'

            for (IX = 1; IX <= NINC; IX++) { // 120
               INCX = INC( IX )
               LX = ABS( INCX )*N

               // Generate the vector X.

               TRANSL = HALF
               smake('GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ), 0, N - 1, RESET, TRANSL );
               if ( N.GT.1 ) {
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               }

               for (IY = 1; IY <= NINC; IY++) { // 110
                  INCY = INC( IY )
                  LY = ABS( INCY )*N

                  // Generate the vector Y.

                  TRANSL = ZERO
                  smake('GE', ' ', ' ', 1, N, Y, 1, YY, ABS( INCY ), 0, N - 1, RESET, TRANSL );
                  if ( N.GT.1 ) {
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  }

                  for (IA = 1; IA <= NALF; IA++) { // 100
                     ALPHA = ALF( IA )
                     NULL = N.LE.0.OR.ALPHA == ZERO

                     // Generate the matrix A.

                     TRANSL = ZERO
                     smake(SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA, LDA, N - 1, N - 1, RESET, TRANSL );

                     NC = NC + 1

                     // Save every datum before calling the subroutine.

                     UPLOS = UPLO
                     NS = N
                     ALS = ALPHA
                     for (I = 1; I <= LAA; I++) { // 10
                        AS( I ) = AA( I )
                     } // 10
                     LDAS = LDA
                     for (I = 1; I <= LX; I++) { // 20
                        XS( I ) = XX( I )
                     } // 20
                     INCXS = INCX
                     for (I = 1; I <= LY; I++) { // 30
                        YS( I ) = YY( I )
                     } // 30
                     INCYS = INCY

                     // Call the subroutine.

                     if ( FULL ) {
                        if (TRACE) WRITE( NTRA, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA;
                        if (REWI) REWIND NTRA;
                        ssyr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA, LDA );
                     } else if ( PACKED ) {
                        if (TRACE) WRITE( NTRA, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY;
                        if (REWI) REWIND NTRA;
                        sspr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA );
                     }

                     // Check if error-exit was taken incorrectly.

                     if ( .NOT.OK ) {
                        WRITE( NOUT, FMT = 9992 )
                        FATAL = true;
                        GO TO 160
                     }

                     // See what data changed inside subroutines.

                     ISAME( 1 ) = UPLO == UPLOS
                     ISAME( 2 ) = NS == N
                     ISAME( 3 ) = ALS == ALPHA
                     ISAME( 4 ) = LSE( XS, XX, LX )
                     ISAME( 5 ) = INCXS == INCX
                     ISAME( 6 ) = LSE( YS, YY, LY )
                     ISAME( 7 ) = INCYS == INCY
                     if ( NULL ) {
                        ISAME( 8 ) = LSE( AS, AA, LAA )
                     } else {
                        ISAME( 8 ) = LSERES( SNAME( 2: 3 ), UPLO, N, N, AS, AA, LDA )
                     }
                     if ( .NOT.PACKED ) {
                        ISAME( 9 ) = LDAS == LDA
                     }

                     // If data was incorrectly changed, report and return.

                     SAME = true;
                     for (I = 1; I <= NARGS; I++) { // 40
                        SAME = SAME.AND.ISAME( I )
                        IF( .NOT.ISAME( I ) ) WRITE( NOUT, FMT = 9998 )I
                     } // 40
                     if ( .NOT.SAME ) {
                        FATAL = true;
                        GO TO 160
                     }

                     if ( .NOT.NULL ) {

                        // Check the result column by column.

                        if ( INCX.GT.0 ) {
                           for (I = 1; I <= N; I++) { // 50
                              Z( I, 1 ) = X( I )
                           } // 50
                        } else {
                           for (I = 1; I <= N; I++) { // 60
                              Z( I, 1 ) = X( N - I + 1 )
                           } // 60
                        }
                        if ( INCY.GT.0 ) {
                           for (I = 1; I <= N; I++) { // 70
                              Z( I, 2 ) = Y( I )
                           } // 70
                        } else {
                           for (I = 1; I <= N; I++) { // 80
                              Z( I, 2 ) = Y( N - I + 1 )
                           } // 80
                        }
                        JA = 1
                        for (J = 1; J <= N; J++) { // 90
                           W( 1 ) = Z( J, 2 )
                           W( 2 ) = Z( J, 1 )
                           if ( UPPER ) {
                              JJ = 1
                              LJ = J
                           } else {
                              JJ = J
                              LJ = N - J + 1
                           }
                           smvch('N', LJ, 2, ALPHA, Z( JJ, 1 ), NMAX, W, 1, ONE, A( JJ, J ), 1, YT, G, AA( JA ), EPS, ERR, FATAL, NOUT, true );
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
                           if (FATAL) GO TO 150;
                        } // 90
                     } else {
                        // Avoid repeating tests with N.le.0.
                        if (N.LE.0) GO TO 140;
                     }

                  } // 100

               } // 110

            } // 120

         } // 130

      } // 140

      // Report result.

      if ( ERRMAX.LT.THRESH ) {
         WRITE( NOUT, FMT = 9999 )SNAME, NC
      } else {
         WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
      }
      GO TO 170

      } // 150
      WRITE( NOUT, FMT = 9995 )J

      } // 160
      WRITE( NOUT, FMT = 9996 )SNAME
      if ( FULL ) {
         WRITE( NOUT, FMT = 9993 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA
      } else if ( PACKED ) {
         WRITE( NOUT, FMT = 9994 )NC, SNAME, UPLO, N, ALPHA, INCX, INCY
      }

      } // 170
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL', 'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH', 'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C', 'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2, ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', Y,', I2, ', AP)                     .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,', I2, ', Y,', I2, ', A,', I3, ')                  .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *', '******' )

      // End of SCHK6

      }
      SUBROUTINE SCHKE( ISNUM, SRNAMT, NOUT )

*  Tests the error exits from the Level 2 Blas.
*  Requires a special version of the error-handling routine XERBLA.
*  ALPHA, BETA, A, X and Y should not need to be defined.

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
      REAL               ALPHA, BETA
      // .. Local Arrays ..
      REAL               A( 1, 1 ), X( 1 ), Y( 1 )
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGBMV, SGEMV, SGER, SSBMV, SSPMV, SSPR, SSPR2, SSYMV, SSYR, SSYR2, STBMV, STBSV, STPMV, STPSV, STRMV, STRSV
      // .. Common blocks ..
      // COMMON /INFOC/INFOT, NOUTC, OK, LERR
      // .. Executable Statements ..
      // OK is set to false by the special version of XERBLA or by CHKXER
      // if anything is wrong.
      OK = true;
      // LERR is set to true by the special version of XERBLA each time
      // it is called, and is then tested and re-set by CHKXER.
      LERR = false;
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160 )ISNUM
   10 INFOT = 1
      sgemv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgemv('N', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      sgemv('N', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      sgemv('N', 2, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      sgemv('N', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      sgemv('N', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   20 INFOT = 1
      sgbmv('/', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgbmv('N', -1, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      sgbmv('N', 0, -1, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgbmv('N', 0, 0, -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      sgbmv('N', 2, 0, 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      sgbmv('N', 0, 0, 1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      sgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 13
      sgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   30 INFOT = 1
      ssymv('/', 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ssymv('U', -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      ssymv('U', 2, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ssymv('U', 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 10
      ssymv('U', 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   40 INFOT = 1
      ssbmv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ssbmv('U', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      ssbmv('U', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      ssbmv('U', 0, 1, ALPHA, A, 1, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      ssbmv('U', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 11
      ssbmv('U', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   50 INFOT = 1
      sspmv('/', 0, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sspmv('U', -1, ALPHA, A, X, 1, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      sspmv('U', 0, ALPHA, A, X, 0, BETA, Y, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      sspmv('U', 0, ALPHA, A, X, 1, BETA, Y, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   60 INFOT = 1
      strmv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      strmv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      strmv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      strmv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      strmv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      strmv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   70 INFOT = 1
      stbmv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      stbmv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      stbmv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      stbmv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      stbmv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      stbmv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      stbmv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   80 INFOT = 1
      stpmv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      stpmv('U', '/', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      stpmv('U', 'N', '/', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      stpmv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      stpmv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
   90 INFOT = 1
      strsv('/', 'N', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      strsv('U', '/', 'N', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      strsv('U', 'N', '/', 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      strsv('U', 'N', 'N', -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 6
      strsv('U', 'N', 'N', 2, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 8
      strsv('U', 'N', 'N', 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  100 INFOT = 1
      stbsv('/', 'N', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      stbsv('U', '/', 'N', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      stbsv('U', 'N', '/', 0, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      stbsv('U', 'N', 'N', -1, 0, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      stbsv('U', 'N', 'N', 0, -1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      stbsv('U', 'N', 'N', 0, 1, A, 1, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      stbsv('U', 'N', 'N', 0, 0, A, 1, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  110 INFOT = 1
      stpsv('/', 'N', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      stpsv('U', '/', 'N', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 3
      stpsv('U', 'N', '/', 0, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 4
      stpsv('U', 'N', 'N', -1, A, X, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      stpsv('U', 'N', 'N', 0, A, X, 0 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  120 INFOT = 1
      sger(-1, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sger(0, -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      sger(0, 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      sger(0, 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      sger(2, 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  130 INFOT = 1
      ssyr('/', 0, ALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ssyr('U', -1, ALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      ssyr('U', 0, ALPHA, X, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ssyr('U', 2, ALPHA, X, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  140 INFOT = 1
      sspr('/', 0, ALPHA, X, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sspr('U', -1, ALPHA, X, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      sspr('U', 0, ALPHA, X, 0, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  150 INFOT = 1
      ssyr2('/', 0, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      ssyr2('U', -1, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      ssyr2('U', 0, ALPHA, X, 0, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      ssyr2('U', 0, ALPHA, X, 1, Y, 0, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 9
      ssyr2('U', 2, ALPHA, X, 1, Y, 1, A, 1 );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      GO TO 170
  160 INFOT = 1
      sspr2('/', 0, ALPHA, X, 1, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 2
      sspr2('U', -1, ALPHA, X, 1, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 5
      sspr2('U', 0, ALPHA, X, 0, Y, 1, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );
      INFOT = 7
      sspr2('U', 0, ALPHA, X, 1, Y, 0, A );
      chkxer(SRNAMT, INFOT, NOUT, LERR, OK );

  170 IF( OK )THEN
         WRITE( NOUT, FMT = 9999 )SRNAMT
      } else {
         WRITE( NOUT, FMT = 9998 )SRNAMT
      }
      RETURN

 9999 FORMAT( ' ', A6, ' PASSED THE TESTS OF ERROR-EXITS' )
 9998 FORMAT( ' ******* ', A6, ' FAILED THE TESTS OF ERROR-EXITS *****', '**' )

      // End of SCHKE

      }
      SUBROUTINE SMAKE( TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, KL, KU, RESET, TRANSL )

*  Generates values for an M by N matrix A within the bandwidth
*  defined by KL and KU.
*  Stores the values in the array AA in the data structure required
*  by the routine, with unwanted elements set to rogue value.

*  TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      REAL               ROGUE
      const              ROGUE = -1.0E10 ;
      // .. Scalar Arguments ..
      REAL               TRANSL
      int                KL, KU, LDA, M, N, NMAX;
      bool               RESET;
      String             DIAG, UPLO;
      String             TYPE;
      // .. Array Arguments ..
      REAL               A( NMAX, * ), AA( * )
      // .. Local Scalars ..
      int                I, I1, I2, I3, IBEG, IEND, IOFF, J, KK;
      bool               GEN, LOWER, SYM, TRI, UNIT, UPPER;
      // .. External Functions ..
      REAL               SBEG
      // EXTERNAL SBEG
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // .. Executable Statements ..
      GEN = TYPE( 1: 1 ) == 'G'
      SYM = TYPE( 1: 1 ) == 'S'
      TRI = TYPE( 1: 1 ) == 'T'
      UPPER = ( SYM.OR.TRI ).AND.UPLO == 'U'
      LOWER = ( SYM.OR.TRI ).AND.UPLO == 'L'
      UNIT = TRI.AND.DIAG == 'U'

      // Generate data in array A.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) ) THEN                IF( ( I.LE.J.AND.J - I.LE.KU ).OR. ( I.GE.J.AND.I - J.LE.KL ) ) {
                  A( I, J ) = SBEG( RESET ) + TRANSL
               } else {
                  A( I, J ) = ZERO
               }
               if ( I != J ) {
                  if ( SYM ) {
                     A( J, I ) = A( I, J )
                  } else if ( TRI ) {
                     A( J, I ) = ZERO
                  }
               }
            }
         } // 10
         if (TRI) A( J, J ) = A( J, J ) + ONE          IF( UNIT ) A( J, J ) = ONE;
      } // 20

      // Store elements in array AS in data structure required by routine.

      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 50
            for (I = 1; I <= M; I++) { // 30
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 30
            for (I = M + 1; I <= LDA; I++) { // 40
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 40
         } // 50
      } else if ( TYPE == 'GB' ) {
         for (J = 1; J <= N; J++) { // 90
            for (I1 = 1; I1 <= KU + 1 - J; I1++) { // 60
               AA( I1 + ( J - 1 )*LDA ) = ROGUE
            } // 60
            DO 70 I2 = I1, MIN( KL + KU + 1, KU + 1 + M - J )
               AA( I2 + ( J - 1 )*LDA ) = A( I2 + J - KU - 1, J )
            } // 70
            for (I3 = I2; I3 <= LDA; I3++) { // 80
               AA( I3 + ( J - 1 )*LDA ) = ROGUE
            } // 80
         } // 90
      } else if ( TYPE == 'SY'.OR.TYPE == 'TR' ) {
         for (J = 1; J <= N; J++) { // 130
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
            for (I = 1; I <= IBEG - 1; I++) { // 100
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 100
            for (I = IBEG; I <= IEND; I++) { // 110
               AA( I + ( J - 1 )*LDA ) = A( I, J )
            } // 110
            for (I = IEND + 1; I <= LDA; I++) { // 120
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 120
         } // 130
      } else if ( TYPE == 'SB'.OR.TYPE == 'TB' ) {
         for (J = 1; J <= N; J++) { // 170
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
            for (I = 1; I <= IBEG - 1; I++) { // 140
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 140
            for (I = IBEG; I <= IEND; I++) { // 150
               AA( I + ( J - 1 )*LDA ) = A( I + J - KK, J )
            } // 150
            for (I = IEND + 1; I <= LDA; I++) { // 160
               AA( I + ( J - 1 )*LDA ) = ROGUE
            } // 160
         } // 170
      } else if ( TYPE == 'SP'.OR.TYPE == 'TP' ) {
         IOFF = 0
         for (J = 1; J <= N; J++) { // 190
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            for (I = IBEG; I <= IEND; I++) { // 180
               IOFF = IOFF + 1
               AA( IOFF ) = A( I, J )
               if ( I == J ) {
                  if (UNIT) AA( IOFF ) = ROGUE;
               }
            } // 180
         } // 190
      }
      RETURN

      // End of SMAKE

      }
      SUBROUTINE SMVCH( TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT, G, YY, EPS, ERR, FATAL, NOUT, MV )

*  Checks the results of the computational tests.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // .. Scalar Arguments ..
      REAL               ALPHA, BETA, EPS, ERR
      int                INCX, INCY, M, N, NMAX, NOUT;
      bool               FATAL, MV;
      String             TRANS;
      // .. Array Arguments ..
      REAL               A( NMAX, * ), G( * ), X( * ), Y( * ), YT( * ), YY( * )
      // .. Local Scalars ..
      REAL               ERRI
      int                I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
      bool               TRAN;
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // .. Executable Statements ..
      TRAN = TRANS == 'T'.OR.TRANS == 'C'
      if ( TRAN ) {
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
      for (I = 1; I <= ML; I++) { // 30
         YT( IY ) = ZERO
         G( IY ) = ZERO
         JX = KX
         if ( TRAN ) {
            for (J = 1; J <= NL; J++) { // 10
               YT( IY ) = YT( IY ) + A( J, I )*X( JX )
               G( IY ) = G( IY ) + ABS( A( J, I )*X( JX ) )
               JX = JX + INCXL
            } // 10
         } else {
            for (J = 1; J <= NL; J++) { // 20
               YT( IY ) = YT( IY ) + A( I, J )*X( JX )
               G( IY ) = G( IY ) + ABS( A( I, J )*X( JX ) )
               JX = JX + INCXL
            } // 20
         }
         YT( IY ) = ALPHA*YT( IY ) + BETA*Y( IY )
         G( IY ) = ABS( ALPHA )*G( IY ) + ABS( BETA*Y( IY ) )
         IY = IY + INCYL
      } // 30

      // Compute the error ratio for this result.

      ERR = ZERO
      for (I = 1; I <= ML; I++) { // 40
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*ABS( INCY ) ) )/EPS
         IF( G( I ) != ZERO ) ERRI = ERRI/G( I )
         ERR = MAX( ERR, ERRI )
         IF( ERR*SQRT( EPS ).GE.ONE ) GO TO 50
      } // 40
      // If the loop completes, all results are at least half accurate.
      GO TO 70

      // Report fatal error.

   50 FATAL = true;
      WRITE( NOUT, FMT = 9999 )
      for (I = 1; I <= ML; I++) { // 60
         if ( MV ) {
            WRITE( NOUT, FMT = 9998 )I, YT( I ), YY( 1 + ( I - 1 )*ABS( INCY ) )
         } else {
            WRITE( NOUT, FMT = 9998 )I, YY( 1 + ( I - 1 )*ABS( INCY ) ), YT(I)
         }
      } // 60

      } // 70
      RETURN

 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL', 'F ACCURATE *******', /'           EXPECTED RESULT   COMPU', 'TED RESULT' )
 9998 FORMAT( 1X, I7, 2G18.6 )

      // End of SMVCH

      }
      bool    FUNCTION LSE( RI, RJ, LR );

*  Tests if two arrays are identical.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LR;
      // .. Array Arguments ..
      REAL               RI( * ), RJ( * )
      // .. Local Scalars ..
      int                I;
      // .. Executable Statements ..
      for (I = 1; I <= LR; I++) { // 10
         IF( RI( I ) != RJ( I ) ) GO TO 20
      } // 10
      LSE = true;
      GO TO 30
      } // 20
      LSE = false;
   30 RETURN

      // End of LSE

      }
      bool    FUNCTION LSERES( TYPE, UPLO, M, N, AA, AS, LDA );

*  Tests if selected elements in two arrays are equal.

*  TYPE is 'GE', 'SY' or 'SP'.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      int                LDA, M, N;
      String             UPLO;
      String             TYPE;
      // .. Array Arguments ..
      REAL               AA( LDA, * ), AS( LDA, * )
      // .. Local Scalars ..
      int                I, IBEG, IEND, J;
      bool               UPPER;
      // .. Executable Statements ..
      UPPER = UPLO == 'U'
      if ( TYPE == 'GE' ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = M + 1; I <= LDA; I++) { // 10
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70
            } // 10
         } // 20
      } else if ( TYPE == 'SY' ) {
         for (J = 1; J <= N; J++) { // 50
            if ( UPPER ) {
               IBEG = 1
               IEND = J
            } else {
               IBEG = J
               IEND = N
            }
            for (I = 1; I <= IBEG - 1; I++) { // 30
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70
            } // 30
            for (I = IEND + 1; I <= LDA; I++) { // 40
               IF( AA( I, J ) != AS( I, J ) ) GO TO 70
            } // 40
         } // 50
      }

      LSERES = true;
      GO TO 80
      } // 70
      LSERES = false;
   80 RETURN

      // End of LSERES

      }
      REAL FUNCTION SBEG( RESET )

*  Generates random numbers uniformly distributed between -0.5 and 0.5.

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.
      // Jeremy Du Croz, NAG Central Office.

      // .. Scalar Arguments ..
      bool               RESET;
      // .. Local Scalars ..
      int                I, IC, MI;
      // .. Save statement ..
      SAVE               I, IC, MI
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // .. Executable Statements ..
      if ( RESET ) {
         // Initialize local variables.
         MI = 891
         I = 7
         IC = 0
         RESET = false;
      }

      // The sequence of values of I is bounded between 1 and 999.
      // If initial I = 1,2,3,6,7 or 9, the period will be 50.
      // If initial I = 4 or 8, the period will be 25.
      // If initial I = 5, the period will be 10.
      // IC is used to break up the period by skipping 1 value of I in 6.

      IC = IC + 1
   10 I = I*MI
      I = I - 1000*( I/1000 )
      if ( IC.GE.5 ) {
         IC = 0
         GO TO 10
      }
      SBEG = REAL( I - 500 )/1001.0
      RETURN

      // End of SBEG

      }
      REAL FUNCTION SDIFF( X, Y )

*  Auxiliary routine for test program for Level 2 Blas.

*  -- Written on 10-August-1987.
      // Richard Hanson, Sandia National Labs.

      // .. Scalar Arguments ..
      REAL               X, Y
      // .. Executable Statements ..
      SDIFF = X - Y
      RETURN

      // End of SDIFF

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
         OK = false;
      }
      LERR = false;
      RETURN

 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ', I2, ' NOT D', 'ETECTED BY ', A6, ' *****' )

      // End of CHKXER

      }
      SUBROUTINE SREGR1( TRANS, M, N, LY, KL, KU, ALPHA, A, LDA, X, INCX, BETA, Y, INCY, YS )

*  Input initialization for regression test.

      // .. Scalar Arguments ..
      String             TRANS;
      int                LY, M, N, KL, KU, LDA, INCX, INCY;
      REAL               ALPHA, BETA
      // .. Array Arguments ..
      REAL               A(LDA,*), X(*), Y(*), YS(*)
      // .. Local Scalars ..
      int                I;
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // .. Executable Statements ..
      TRANS = 'T'
      M = 0
      N = 5
      KL = 0
      KU = 0
      ALPHA = 1.0
      LDA = MAX( 1, M )
      INCX = 1
      BETA = -0.7
      INCY = 1
      LY = ABS( INCY )*N
      for (I = 1; I <= LY; I++) { // 10
         Y( I ) = 42.0 + REAL( I )
         YS( I ) = Y( I )
      } // 10
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
      // COMMON /INFOC/INFOT, NOUT, OK, LERR
      // COMMON /SRNAMC/SRNAMT
      // .. Executable Statements ..
      LERR = true;
      if ( INFO != INFOT ) {
         if ( INFOT != 0 ) {
            WRITE( NOUT, FMT = 9999 )INFO, INFOT
         } else {
            WRITE( NOUT, FMT = 9997 )INFO
         }
         OK = false;
      }
      if ( SRNAME != SRNAMT ) {
         WRITE( NOUT, FMT = 9998 )SRNAME, SRNAMT
         OK = false;
      }
      RETURN

 9999 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' INSTEAD', ' OF ', I2, ' *******' )
 9998 FORMAT( ' ******* XERBLA WAS CALLED WITH SRNAME = ', A6, ' INSTE', 'AD OF ', A6, ' *******' )
 9997 FORMAT( ' ******* XERBLA WAS CALLED WITH INFO = ', I6, ' *******' )

      // End of XERBLA

      }
