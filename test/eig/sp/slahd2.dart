      void slahd2(IOUNIT, PATH ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                IOUNIT;
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               CORZ, SORD;
      String             C2;
      int                J;
      // ..
      // .. External Functions ..
      //- bool               lsame, LSAMEN;
      // EXTERNAL lsame, LSAMEN

      if (IOUNIT <= 0) return;
      SORD = lsame( PATH, 'S' ) || lsame( PATH, 'D' );
      CORZ = lsame( PATH, 'C' ) || lsame( PATH, 'Z' );
      if ( !SORD && !CORZ ) {
         WRITE( IOUNIT, FMT = 9999 )PATH;
      }
      C2 = PATH( 2: 3 );

      if ( lsamen( 2, C2, 'HS' ) ) {
         if ( SORD ) {

            // Real Non-symmetric Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9998 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9988 );
            WRITE( IOUNIT, FMT = 9987 );
            WRITE( IOUNIT, FMT = 9986 )'pairs ', 'pairs ', 'prs.', 'prs.';
            WRITE( IOUNIT, FMT = 9985 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9984 )'orthogonal', '''=transpose', ( '''', J = 1, 6 );

         } else {

            // Complex Non-symmetric Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9997 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9988 );
            WRITE( IOUNIT, FMT = 9987 );
            WRITE( IOUNIT, FMT = 9986 )'e.vals', 'e.vals', 'e.vs', 'e.vs';
            WRITE( IOUNIT, FMT = 9985 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9984 )'unitary', '*=conj.transp.', ( '*', J = 1, 6 );
         }

      } else if ( lsamen( 2, C2, 'ST' ) ) {

         if ( SORD ) {

            // Real Symmetric Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9996 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9983 );
            WRITE( IOUNIT, FMT = 9982 );
            WRITE( IOUNIT, FMT = 9981 )'Symmetric';

            // Tests performed

            WRITE( IOUNIT, FMT = 9968 );

         } else {

            // Complex Hermitian Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9995 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9983 );
            WRITE( IOUNIT, FMT = 9982 );
            WRITE( IOUNIT, FMT = 9981 )'Hermitian';

            // Tests performed

            WRITE( IOUNIT, FMT = 9967 );
         }

      } else if ( lsamen( 2, C2, 'SG' ) ) {

         if ( SORD ) {

            // Real Symmetric Generalized Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9992 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9980 );
            WRITE( IOUNIT, FMT = 9979 );
            WRITE( IOUNIT, FMT = 9978 )'Symmetric';

            // Tests performed

            WRITE( IOUNIT, FMT = 9977 );
            WRITE( IOUNIT, FMT = 9976 );

         } else {

            // Complex Hermitian Generalized Eigenvalue Problem:

            WRITE( IOUNIT, FMT = 9991 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9980 );
            WRITE( IOUNIT, FMT = 9979 );
            WRITE( IOUNIT, FMT = 9978 )'Hermitian';

            // Tests performed

            WRITE( IOUNIT, FMT = 9975 );
            WRITE( IOUNIT, FMT = 9974 );

         }

      } else if ( lsamen( 2, C2, 'BD' ) ) {

         if ( SORD ) {

            // Real Singular Value Decomposition:

            WRITE( IOUNIT, FMT = 9994 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9973 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9972 )'orthogonal';
            WRITE( IOUNIT, FMT = 9971 );
         } else {

            // Complex Singular Value Decomposition:

            WRITE( IOUNIT, FMT = 9993 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9973 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9972 )'unitary   ';
            WRITE( IOUNIT, FMT = 9971 );
         }

      } else if ( lsamen( 2, C2, 'BB' ) ) {

         if ( SORD ) {

            // Real General Band reduction to bidiagonal form:

            WRITE( IOUNIT, FMT = 9990 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9970 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9969 )'orthogonal';
         } else {

            // Complex Band reduction to bidiagonal form:

            WRITE( IOUNIT, FMT = 9989 )PATH;

            // Matrix types

            WRITE( IOUNIT, FMT = 9970 );

            // Tests performed

            WRITE( IOUNIT, FMT = 9969 )'unitary   ';
         }

      } else {

         WRITE( IOUNIT, FMT = 9999 )PATH;
         return;
      }

      return;

 9999 FORMAT(' ${.a3}:  no header available' );
 9998 FORMAT( / 1X, '${.a3} -- Real Non-symmetric eigenvalue problem' );
 9997 FORMAT( / 1X, '${.a3} -- Complex Non-symmetric eigenvalue problem' );
 9996 FORMAT( / 1X, '${.a3} -- Real Symmetric eigenvalue problem' );
 9995 FORMAT( / 1X, '${.a3} -- Complex Hermitian eigenvalue problem' );
 9994 FORMAT( / 1X, '${.a3} -- Real Singular Value Decomposition' );
 9993 FORMAT( / 1X, '${.a3} -- Complex Singular Value Decomposition' );
 9992 FORMAT( / 1X, '${.a3} -- Real Symmetric Generalized eigenvalue problem' );
 9991 FORMAT( / 1X, '${.a3} -- Complex Hermitian Generalized eigenvalue problem' );
 9990 FORMAT( / 1X, '${.a3} -- Real Band reduc. to bidiagonal form' );
 9989 FORMAT( / 1X, '${.a3} -- Complex Band reduc. to bidiagonal form' );

 9988 FORMAT( ' Matrix types (see xCHKHS for details): ' );

 9987 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.' );
 9986 FORMAT( ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex ', A6, / ' 12=Well-cond., random complex ${.a6}    17=Ill-cond., large rand. complx ', A4, / ' 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ', A4 );
 9985 FORMAT( ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.   ' );
 9984 FORMAT('\n Tests performed:   (H is Hessenberg, T is Schur, U and Z are ${},', / 20X, A, ', W is a diagonal matrix of eigenvalues,\n${' ' * 20}L and R are the left and right eigenvector matrices)\n  1 = | A - U H U${.a1} | / ( |A| n ulp )           2 = | I - U U${.a1} | / ( n ulp )\n  3 = | H - Z T Z${.a1} | / ( |H| n ulp )           4 = | I - Z Z${.a1} | / ( n ulp )\n  5 = | A - UZ T (UZ)${.a1} | / ( |A| n ulp )       6 = | I - UZ (UZ)${.a1} | / ( n ulp )\n  7 = | T(e.vects.) - T(no e.vects.) | / ( |T| ulp )\n  8 = | W(e.vects.) - W(no e.vects.) | / ( |W| ulp )\n  9 = | TR - RW | / ( |T| |R| ulp )      10 = | LT - WL | / ( |T| |L| ulp )\n 11= |HX - XW| / (|H| |X| ulp)  (inv.it) 12= |YH - WY| / (|H| |Y| ulp)  (inv.it)' );

      // Symmetric/Hermitian eigenproblem

 9983 FORMAT( ' Matrix types (see xDRVST for details): ' );

 9982 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.' );
 9981 FORMAT( ' Dense ${} Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.' );

      // Symmetric/Hermitian Generalized eigenproblem

 9980 FORMAT( ' Matrix types (see xDRVSG for details): ' );

 9979 FORMAT('\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.' );
 9978 FORMAT( ' Dense or Banded ${} Matrices: \n  8=Evenly spaced eigenvals.          15=Matrix with small random entries.\n  9=Geometrically spaced eigenvals.   16=Evenly spaced eigenvals, KA=1, KB=1.\n 10=Clustered eigenvalues.            17=Evenly spaced eigenvals, KA=2, KB=1.\n 11=Large, evenly spaced eigenvals.   18=Evenly spaced eigenvals, KA=2, KB=2.\n 12=Small, evenly spaced eigenvals.   19=Evenly spaced eigenvals, KA=3, KB=1.\n 13=Matrix with random O(1) entries.  20=Evenly spaced eigenvals, KA=3, KB=2.\n 14=Matrix with large random entries. 21=Evenly spaced eigenvals, KA=3, KB=3.' );
 9977 FORMAT('\n Tests performed:   \n( For each pair (A,B), where A is of the given type \n and B is a random well-conditioned matrix. D is \n diagonal, and Z is orthogonal. )\n 1 = SSYGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 2 = SSPGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 3 = SSBGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 4 = SSYGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 5 = SSPGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 6 = SSBGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     ' );
 9976 FORMAT( ' 7 = SSYGV, with ITYPE=2 and UPLO=''U'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 8 = SSPGV, with ITYPE=2 and UPLO=''U'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 9 = SSPGV, with ITYPE=2 and UPLO=''L'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n10 = SSPGV, with ITYPE=2 and UPLO=''L'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n11 = SSYGV, with ITYPE=3 and UPLO=''U'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n12 = SSPGV, with ITYPE=3 and UPLO=''U'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n13 = SSYGV, with ITYPE=3 and UPLO=''L'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n14 = SSPGV, with ITYPE=3 and UPLO=''L'':  | B A Z - Z D | / ( |A| |Z| n ulp )     ' );
 9975 FORMAT('\n Tests performed:   \n( For each pair (A,B), where A is of the given type \n and B is a random well-conditioned matrix. D is \n diagonal, and Z is unitary. )\n 1 = CHEGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 2 = CHPGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 3 = CHBGV, with ITYPE=1 and UPLO=''U'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 4 = CHEGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 5 = CHPGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     \n 6 = CHBGV, with ITYPE=1 and UPLO=''L'':  | A Z - B Z D | / ( |A| |Z| n ulp )     ' );
 9974 FORMAT( ' 7 = CHEGV, with ITYPE=2 and UPLO=''U'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 8 = CHPGV, with ITYPE=2 and UPLO=''U'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n 9 = CHPGV, with ITYPE=2 and UPLO=''L'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n10 = CHPGV, with ITYPE=2 and UPLO=''L'':  | A B Z - Z D | / ( |A| |Z| n ulp )     \n11 = CHEGV, with ITYPE=3 and UPLO=''U'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n12 = CHPGV, with ITYPE=3 and UPLO=''U'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n13 = CHEGV, with ITYPE=3 and UPLO=''L'':  | B A Z - Z D | / ( |A| |Z| n ulp )     \n14 = CHPGV, with ITYPE=3 and UPLO=''L'':  | B A Z - Z D | / ( |A| |Z| n ulp )     ' );

      // Singular Value Decomposition

 9973 FORMAT( ' Matrix types (see xCHKBD for details):\n Diagonal matrices:\n   1: Zero${' ' * 28} 5: Clustered entries\n   2: Identity${' ' * 24} 6: Large, evenly spaced entries\n   3: Evenly spaced entries${' ' * 11} 7: Small, evenly spaced entries\n   4: Geometrically spaced entries\n General matrices:\n   8: Evenly spaced sing. vals.${' ' * 7}12: Small, evenly spaced sing vals\n   9: Geometrically spaced sing vals  13: Random, O(1) entries\n  10: Clustered sing. vals.${' ' * 11}14: Random, scaled near overflow\n  11: Large, evenly spaced sing vals  15: Random, scaled near underflow' );

 9972 FORMAT('\n Test ratios:  (B: bidiagonal, S: diagonal, Q, P, U, and V: ', A10, / 16X, 'X: m x nrhs, Y = Q'' X, and Z = U'' Y)' );
 9971 FORMAT( '   1: norm( A - Q B P'' ) / ( norm(A) max(m,n) ulp )\n   2: norm( I - Q'' Q )   / ( m ulp )\n   3: norm( I - P'' P )   / ( n ulp )\n   4: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )\n   5: norm( Y - U Z )    / ( norm(Z) max(min(m,n),k) ulp )\n   6: norm( I - U'' U )   / ( min(m,n) ulp )\n   7: norm( I - V'' V )   / ( min(m,n) ulp )\n   8: Test ordering of S  (0 if nondecreasing, 1/ulp  otherwise)\n   9: norm( S - S1 )     / ( norm(S) ulp ), where S1 is computed\n${' ' * 43} without computing U and V''\n  10: Sturm sequence test (0 if sing. vals of B within THRESH of S)\n  11: norm( A - (QU) S (V'' P'') ) / ( norm(A) max(m,n) ulp )\n  12: norm( X - (QU) Z )         / ( |X| max(M,k) ulp )\n  13: norm( I - (QU)''(QU) )      / ( M ulp )\n  14: norm( I - (V'' P'') (P V) )  / ( N ulp )\n  15: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )\n  16: norm( I - U'' U )   / ( min(m,n) ulp )\n  17: norm( I - V'' V )   / ( min(m,n) ulp )',/ '  18: Test ordering of S  (0 if nondecreasing, 1/ulp ',  ' otherwise)',/ '  19: norm( S - S1 )     / ( norm(S) ulp ),',  ' where S1 is computed', / 43X,  ' without computing U and V''',/ '  20: norm( B - U S V'' )  / ( norm(B) min(m,n) ulp )',  '  SBDSVX(V,A)',/ '  21: norm( I - U'' U )    / ( min(m,n) ulp )',/ '  22: norm( I - V'' V )    / ( min(m,n) ulp )',/ '  23: Test ordering of S  (0 if nondecreasing, 1/ulp ',  ' otherwise)',/ '  24: norm( S - S1 )      / ( norm(S) ulp ),',  ' where S1 is computed', / 44X,  ' without computing U and V''',/ '  25: norm( S - U'' B V ) / ( norm(B) n ulp )',  '  SBDSVX(V,I)',/ '  26: norm( I - U'' U )    / ( min(m,n) ulp )',/ '  27: norm( I - V'' V )    / ( min(m,n) ulp )',/ '  28: Test ordering of S  (0 if nondecreasing, 1/ulp ',  ' otherwise)',/ '  29: norm( S - S1 )      / ( norm(S) ulp ),',  ' where S1 is computed', / 44X,  ' without computing U and V''',/ '  30: norm( S - U'' B V ) / ( norm(B) n ulp )',  '  SBDSVX(V,V)',/ '  31: norm( I - U'' U )    / ( min(m,n) ulp )',/ '  32: norm( I - V'' V )    / ( min(m,n) ulp )',/ '  33: Test ordering of S  (0 if nondecreasing, 1/ulp ',  ' otherwise)',/ '  34: norm( S - S1 )      / ( norm(S) ulp ),',  ' where S1 is computed', / 44X,  ' without computing U and V''' );

      // Band reduction to bidiagonal form

 9970 FORMAT( ' Matrix types (see xCHKBB for details):\n Diagonal matrices:\n   1: Zero${' ' * 28} 5: Clustered entries\n   2: Identity${' ' * 24} 6: Large, evenly spaced entries\n   3: Evenly spaced entries${' ' * 11} 7: Small, evenly spaced entries\n   4: Geometrically spaced entries\n General matrices:\n   8: Evenly spaced sing. vals.${' ' * 7}12: Small, evenly spaced sing vals\n   9: Geometrically spaced sing vals  13: Random, O(1) entries\n  10: Clustered sing. vals.${' ' * 11}14: Random, scaled near overflow\n  11: Large, evenly spaced sing vals  15: Random, scaled near underflow' );

 9969 FORMAT('\n Test ratios:  (B: upper bidiagonal, Q and P: ', A10, / 16X, 'C: m x nrhs, PT = P'', Y = Q'' C)\n 1: norm( A - Q B PT ) / ( norm(A) max(m,n) ulp )\n 2: norm( I - Q'' Q )   / ( m ulp )\n 3: norm( I - PT PT'' )   / ( n ulp )\n 4: norm( Y - Q'' C )   / ( norm(Y) max(m,nrhs) ulp )' );
 9968 FORMAT('\n Tests performed:  See sdrvst.f' );
 9967 FORMAT('\n Tests performed:  See cdrvst.f' );
      }
