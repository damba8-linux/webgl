//WebGL kontekst
function getWebGL(kanID) {
  var kan = document.getElementById(kanID);
  if (kan == null) {
    alert("Nemate canvas na stranici.");
    return;
  }
  var names = ["webgl", "experimental-webgl", "webkit-3d", "moz-webgl"];
  var ctx = null;
  for (var i = 0; i < names.length; i++) {
    try {
      ctx = kan.getContext(names[i]);
    }
    catch(e){
      if (ctx) break;
    }
    if (ctx == null) alert("WebGL nije dostupan.");
    else return ctx;
  }
}

//kreiranje shadera
function makeShader(gl, source, type) {
  //napravi shader objekt
  var shader = gl.createShader(type);
  //pridruzi shader source kod
  gl.shaderSource(shader, source);
  //kompajliraj shader
  gl.compileShader(shader);
  //provjeri je li sve ok
  var uspjeh = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
  if (!uspjeh) {
    throw "Shader nije kompajliran: " + gl.getShaderInfoLog(shader);
  }
  return shader;
}

//povezivanje shadera u program
function makeProgram(gl, vshader, fshader) {
  var program = gl.createProgram();
  //pridruzi shadere
  gl.attachShader(program, vshader);
  gl.attachShader(program, fshader);
  //povezi shadere u program
  gl.linkProgram(program);
  //provjeri je li dobro povezano
  var uspjeh = gl.getProgramParameter(program, gl.LINK_STATUS);
  if (!uspjeh) {
    throw "Program nije kreiran kako treba: " + gl.getProgramInfoLog(program);
  }
  return program;
}

//kreiranje shadera iz <script> taga
function kompajlirajShader(gl, id, type) {
  //potrazi skriptu u dokumentu
  var shaderSkripta = document.getElementById(id);
  if (!shaderSkripta) {
    throw "Nepoznata skripta: " + id;
  }
  //uzmi sadrzaj skripte
  var shaderSource = shaderSkripta.text.trim();
  //ako nismo specificirali tip shadera
  if (!type) {
    if (shaderSkripta.type == "x-shader/x-vertex") {
      type = gl.VERTEX_SHADER;
    } else if (shaderSkripta.type == "x-shader/x-fragment") {
      type = gl.FRAGMENT_SHADER;
    } else if (!type) {
      throw "Nije specificiran tip shadera."
    }
  }
  return makeShader(gl, shaderSource, type);
}

//napravi sve zajedno
function napraviProgram(gl, vsID, fsID) {
  var vertexShader = kompajlirajShader(gl, vsID);
  var fragmentShader = kompajlirajShader(gl, fsID);
  return makeProgram(gl, vertexShader, fragmentShader);
}
