function MT3D() {
  'use strict';
  var matrica = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
  var kamera = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
  var projekcija = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
  var transformacije = [];
  var pamtiSkaliranje = [];
  var skaliranje = false; //da li je primijenjeno skaliranje; olaksava racunanje matrice za transformaciju normale
  //vracanje  matrice transformacije (za prvih pet labosa)
  this.getMatrica = function() {
    return matrica;
  };
  //pretvori kamera*matrica u polje
  this.poljeModel = function() {
    var rez = this.mnoziMatrice(kamera, matrica);
    var v = [];
    for(var i = 0; i < 4; i++) {
      for (var j = 0; j < 4; j++) {
        v.push(rez[j][i]);
      }
    }
    return v;
  };
  //pretvori matricu kamere u polje
  this.poljeKamera = function() {
    var v = [];
    for(var i = 0; i < 4; i++) {
      for (var j = 0; j < 4; j++) {
        v.push(kamera[j][i]);
      }
    }
    return v;
  };
  //pretvori projekcija*kamera*matrica u polje
  this.poljeMatrica = function() {
    var rez = this.mnoziMatrice(projekcija, this.mnoziMatrice(kamera, matrica));
    var v = [];
    for(var i = 0; i < 4; i++) {
      for (var j = 0; j < 4; j++) {
        v.push(rez[j][i]);
      }
    }
    return v;
  };
  //pretvori matricu za transformaciju normala u polje
  this.poljeNormala = function() {
    var m = [];
    var matNormala;
    var v = [];
    var rez = this.mnoziMatrice(kamera, matrica);
    //izvuci 3x3 podmatricu
    for (var i = 0; i < 3; i++) {
      m[i] = [];
      for (var j = 0; j < 3; j++) {
        m[i][j] = rez[i][j];
      }
    }
    if (!skaliranje) matNormala = m;
    else {
      var det = 1 / (m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] -
                m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1]);
      var M00 = det*(m[1][1]*m[2][2]-m[2][1]*m[1][2]);
      var M01 = det*(m[2][0]*m[1][2]-m[1][0]*m[2][2]);
      var M02 = det*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
      var M10 = det*(m[0][2]*m[2][1]-m[0][1]*m[2][2]);
      var M11 = det*(m[0][0]*m[2][2]-m[0][2]*m[2][0]);
      var M12 = det*(m[0][1]*m[2][0]-m[0][0]*m[2][1]);
      var M20 = det*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
      var M21 = det*(m[0][2]*m[1][0]-m[0][0]*m[1][2]);
      var M22 = det*(m[0][0]*m[1][1]-m[0][1]*m[1][0]);
      matNormala = [[M00, M01, M02], [M10, M11, M12], [M20, M21, M22]];
    }
    for(var i = 0; i < 3; i++) {
      for (var j = 0; j < 3; j++) {
        v.push(matNormala[j][i]);
      }
    }
    return v;
  };
  //pretvori matricu za transformaciju normala u polje prikladno za Uniform buffer (webgl2)
  this.poljeNormalaBuffer = function() {
    var v = this.poljeNormala();
    v.splice(3,0,0);
    v.splice(7,0,0);
    v.push(0);
    return v; 
  };
  //spremi trenutnu transformaciju na stog
  this.spremiMatricu = function() {
    if (transformacije.length > 32) {
      console.log("Error: stack is full.");
    } else {
      transformacije.push(matrica);
      pamtiSkaliranje.push(skaliranje);
    }
  };
  //vrati zadnje spremljenu transformaciju sa stoga
  this.vratiMatricu = function() {
    if (transformacije.length == 0) {
      console.log("Error: stack is empty.");
    } else {
      matrica = transformacije.pop();
      skaliranje = pamtiSkaliranje.pop();
    }
  };
  //mnozenje matrica
  this.mnoziMatrice = function(m1, m2) {
    var i, j, k;
    var rez = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          rez[i][j] += m1[i][k] * m2[k][j];
        }
      }
    }
    return rez;
  };
  //kompozicija transformacija
  this.mult = function(m) {
    var i, j, k;
    var m1 = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          m1[i][j] += matrica[i][k] * m[k][j];
        }
      }
    }
    matrica = m1;
  };
  //slika tocke
  this.slika_tocke = function(tocka) {
    var ImTocka = [0,0,0,0];
    var i, k;
    for (i = 0; i < 4; i++) {
      for (k = 0; k < 4; k++) {
        ImTocka[i] += matrica[i][k] * tocka[k];
      }
    }
    return ImTocka;
  };
  //slika tocke primjenom kamere
  this.kamera_slika_tocke = function(tocka) {
    var ImTocka = [0,0,0,0];
    var i, k;
    for (i = 0; i < 4; i++) {
      for (k = 0; k < 4; k++) {
        ImTocka[i] += kamera[i][k] * tocka[k];
      }
    }
    return ImTocka;
  };
  //vektorski produkt
  this.VP = function(u, v) {
    var vek = [0,0,0];
    vek[0] = u[1]*v[2] - u[2]*v[1];
    vek[1] = u[2]*v[0] - u[0]*v[2];
    vek[2] = u[0]*v[1] - u[1]*v[0];
    return vek;
  };
  //translacija
  this.pomakni = function(px, py, pz) {
    var m = [[1,0,0,px],[0,1,0,py],[0,0,1,pz],[0,0,0,1]];
    this.mult(m);
  };
  //skaliranje
  this.skaliraj = function(sx, sy, sz) {
    var m = [[sx,0,0,0],[0,sy,0,0],[0,0,sz,0],[0,0,0,1]];
    if (!skaliranje) skaliranje = true;
    this.mult(m);
  };
  //zrcaljenja na koordinatnim osima
  this.zrcaliNaX = function() {
    var m = [[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]];
    this.mult(m);
  };
  this.zrcaliNaY = function() {
    var m = [[-1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]];
    this.mult(m);
  };
  this.zrcaliNaZ = function() {
    var m = [[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]];
    this.mult(m);
  };
  //zrcaljenja na koordinatnim ravninama
  this.zrcaliNaXY = function() {
    var m = [[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]];
    this.mult(m);
  }
  this.zrcaliNaXZ = function() {
    var m = [[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]];
    this.mult(m);
  }
  this.zrcaliNaYZ = function() {
    var m = [[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
    this.mult(m);
  }
  //rotacije oko koordinatnih osi u stupnjevima
  this.rotirajX = function(kut) {
    var fi = kut * Math.PI / 180;
    var cosfi = Math.cos(fi);
    var sinfi = Math.sin(fi);
    var m = [[1,0,0,0],[0,cosfi,-sinfi,0],[0,sinfi,cosfi,0],[0,0,0,1]];
    this.mult(m);
  };
  this.rotirajY = function(kut) {
    var fi = kut * Math.PI / 180;
    var cosfi = Math.cos(fi);
    var sinfi = Math.sin(fi);
    var m = [[cosfi,0,sinfi,0],[0,1,0,0],[-sinfi,0,cosfi,0],[0,0,0,1]];
    this.mult(m);
  };
  this.rotirajZ = function(kut) {
    var fi = kut * Math.PI / 180;
    var cosfi = Math.cos(fi);
    var sinfi = Math.sin(fi);
    var m = [[cosfi,-sinfi,0,0],[sinfi,cosfi,0,0],[0,0,1,0],[0,0,0,1]];
    this.mult(m);
  };
  //identiteta
  this.identitet = function() {
    matrica = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
    skaliranje = false;
  };
  //smicanje
  this.smicanje = function(alpha, beta, gamma) {
    var tanA = Math.tan(alpha*Math.PI/180);
    var tanB = Math.tan(beta*Math.PI/180);
    var tanC = Math.tan(gamma*Math.PI/180);
    var m = [[1,tanB,tanC,0],[tanA,1,tanC,0],[tanA,tanB,1,0],[0,0,0,1]];
    if (!skaliranje) skaliranje = true;
    this.mult(m);
  };
  //opcenitije transformacije
  this.rotiraj_oko_osi = function(x0, y0, z0, u1, u2, u3, kut) {
    // os je zadana tockom (x0,y0,z0) i vektorom smjera (u1,u2,u3)
    var korijen = Math.sqrt(u1*u1+u2*u2+u3*u3);
    var a = u1/korijen;
    var b = u2/korijen;
    var c = u3/korijen;
    var d = Math.sqrt(b*b+c*c);
    var Rx1 = [[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]];
    var Ry1 = [[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]];
    var Rx2 = [[1,0,0,0],[0,c/d,b/d,0],[0,-b/d,c/d,0],[0,0,0,1]];
    var Ry2 = [[d,0,a,0],[0,1,0,0],[-a,0,d,0],[0,0,0,1]];
    this.pomakni(x0,y0,z0);
    this.mult(Rx2);
    this.mult(Ry2);
    this.rotirajZ(kut);
    this.mult(Ry1);
    this.mult(Rx1);
    this.pomakni(-x0,-y0,-z0);
  };
  this.zrcali_os = function(x0, y0, z0, u1, u2, u3) {
    // os je zadana tockom (x0,y0,z0) i vektorom smjera (u1,u2,u3)
    var korijen = Math.sqrt(u1*u1+u2*u2+u3*u3);
    var a = u1/korijen;
    var b = u2/korijen;
    var c = u3/korijen;
    var d = Math.sqrt(b*b+c*c);
    var Rx1 = [[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]];
    var Ry1 = [[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]];
    var Rx2 = [[1,0,0,0],[0,c/d,b/d,0],[0,-b/d,c/d,0],[0,0,0,1]];
    var Ry2 = [[d,0,a,0],[0,1,0,0],[-a,0,d,0],[0,0,0,1]];
    this.pomakni(x0,y0,z0);
    this.mult(Rx2);
    this.mult(Ry2);
    this.zrcaliNaZ();
    this.mult(Ry1);
    this.mult(Rx1);
    this.pomakni(-x0,-y0,-z0);
  };
  this.zrcali_ravnina = function(x0, y0, z0, u1, u2, u3) {
    // ravnina je zadana tockom (x0,y0,z0) i normalom (u1,u2,u3)
    var korijen = Math.sqrt(u1*u1+u2*u2+u3*u3);
    var a = u1/korijen;
    var b = u2/korijen;
    var c = u3/korijen;
    var d = Math.sqrt(b*b+c*c);
    var Rx1 = [[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]];
    var Ry1 = [[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]];
    var Rx2 = [[1,0,0,0],[0,c/d,b/d,0],[0,-b/d,c/d,0],[0,0,0,1]];
    var Ry2 = [[d,0,a,0],[0,1,0,0],[-a,0,d,0],[0,0,0,1]];
    this.pomakni(x0,y0,z0);
    this.mult(Rx2);
    this.mult(Ry2);
    this.zrcaliNaXY();
    this.mult(Ry1);
    this.mult(Rx1);
    this.pomakni(-x0,-y0,-z0);
  };
  this.zrcali_ravnina2 = function(A, B, C, D) {
    // ravnina je zadana u opcem obliku Ax+By+Cz+D=0
    var x0, y0, z0;
    if (A!=0) {
      x0 = -D/A;
      y0 = 0;
      z0 = 0;
    } else if (B!=0) {
      x0 = 0;
      y0 = -D/B;
      z0 = 0;
    } else {
      x0 = 0;
      y0 = 0;
      z0 = -D/C;
    }
    var korijen = Math.sqrt(A*A+B*B+C*C);
    var a = A/korijen;
    var b = B/korijen;
    var c = C/korijen;
    var d = Math.sqrt(b*b+c*c);
    var Rx1 = [[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]];
    var Ry1 = [[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]];
    var Rx2 = [[1,0,0,0],[0,c/d,b/d,0],[0,-b/d,c/d,0],[0,0,0,1]];
    var Ry2 = [[d,0,a,0],[0,1,0,0],[-a,0,d,0],[0,0,0,1]];
    this.pomakni(x0,y0,z0);
    this.mult(Rx2);
    this.mult(Ry2);
    this.zrcaliNaXY();
    this.mult(Ry1);
    this.mult(Rx1);
    this.pomakni(-x0,-y0,-z0);
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //i pritom se ne cuvaju proporcije projekcijom na canvas
  this.OrtogonalnaProjekcija = function(xmin, xmax, ymin, ymax, zmin, zmax) {
    projekcija = [[2/(xmax-xmin),0,0,(xmin+xmax)/(xmin-xmax)],
             [0,2/(ymax-ymin),0,(ymin+ymax)/(ymin-ymax)],
             [0,0,2/(zmin-zmax),(zmin+zmax)/(zmin-zmax)],
             [0,0,0,1]];
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //i pritom se cuvaju proporcije projekcijom na canvas tako da se po potrebi
  //poveca ili smanji interval [ymin,ymax]
  this.OrtogonalnaProjekcijaX = function(xmin, xmax, ymin, ymax, zmin, zmax, w, h) {
    var k = (h / w * (xmax - xmin) - (ymax - ymin)) / 2;
    var y1 = ymin - k;
    var y2 = ymax + k;
    this.OrtogonalnaProjekcija(xmin, xmax, y1, y2, zmin, zmax);
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //i pritom se cuvaju proporcije projekcijom na canvas tako da se po potrebi
  //poveca ili smanji interval [xmin,xmax]
  this.OrtogonalnaProjekcijaY = function(xmin, xmax, ymin, ymax, zmin, zmax, w, h) {
    var k = (w / h * (ymax - ymin) - (xmax - xmin)) / 2;
    var x1 = xmin - k;
    var x2 = xmax + k;
    this.OrtogonalnaProjekcija(x1, x2, ymin, ymax, zmin, zmax);
  };
  //koordinatni sustav kamere
  this.postaviKameru = function(x0, y0, z0, x1, y1, z1, Vx, Vy, Vz) {
    var V = [Vx, Vy, Vz];
    var normaN = Math.sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1));
    var n = [(x0-x1)/normaN, (y0-y1)/normaN, (z0-z1)/normaN];
    var U = this.VP(V, n);
    var normaU = Math.sqrt(U[0]*U[0]+U[1]*U[1]+U[2]*U[2]);
    var u = [U[0]/normaU, U[1]/normaU, U[2]/normaU];
    var v = this.VP(n, u);
    var mtr = [[u[0], u[1], u[2], -u[0]*x0-u[1]*y0-u[2]*z0],
               [v[0], v[1], v[2], -v[0]*x0-v[1]*y0-v[2]*z0],
               [n[0], n[1], n[2], -n[0]*x0-n[1]*y0-n[2]*z0],
               [0, 0, 0, 1] ];
    kamera = mtr;
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //tako da odbacivanjem z-koordinate dobijemo perspektivu
  //ne cuvaju se proprorcije na canvas
  this.PerspektivnaProjekcija = function(xmin, xmax, ymin, ymax, zmin, zmax) {
    projekcija = [[2*zmin/(xmax-xmin),0,(xmax+xmin)/(xmax-xmin),0],
             [0,2*zmin/(ymax-ymin),(ymax+ymin)/(ymax-ymin),0],
             [0,0,(zmin+zmax)/(zmin-zmax),2*zmin*zmax/(zmin-zmax)],
             [0,0,-1,0]];
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //tako da odbacivanjem z-koordinate dobijemo perspektivu
  //cuvaju se proporcije projekcijom na canvas tako da se po potrebi
  //poveca ili smanji interval [ymin,ymax]
  this.PerspektivnaProjekcijaX = function(xmin, xmax, ymin, ymax, zmin, zmax, w, h) {
    var k = (h / w * (xmax - xmin) - (ymax - ymin)) / 2;
    var y1 = ymin - k;
    var y2 = ymax + k;
    this.PerspektivnaProjekcija(xmin, xmax, y1, y2, zmin, zmax);
  };
  //preslikava [xmin,xmax] x [ymin,ymax] x [zmin,zmax] u [-1,1] x [-1,1] x [-1,1]
  //tako da odbacivanjem z-koordinate dobijemo perspektivu
  //cuvaju se proporcije projekcijom na canvas tako da se po potrebi
  //poveca ili smanji interval [xmin,xmax]
  this.PerspektivnaProjekcijaY = function(xmin, xmax, ymin, ymax, zmin, zmax, w, h) {
    var k = (w / h * (ymax - ymin) - (xmax - xmin)) / 2;
    var x1 = xmin - k;
    var x2 = xmax + k;
    this.PerspektivnaProjekcija(x1, x2, ymin, ymax, zmin, zmax);
  };
  //perspektivnu projekciju zadajemo s vertikalnim kutom pogleda,
  //omjerom sirine i visine slike, te granicama na z-osi
  this.Perspektiva = function(theta, omjer, near, far) {
    var top = near * Math.tan(theta * Math.PI / 360 );
    var right = omjer * top;
    this.PerspektivnaProjekcija(-right, right, -top, top, near, far);
  };
}
