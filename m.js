(function(window){
	'use_strict'

	function define_minimal(){
		var m = {};
		var PI = 3.14159265359;

		m.greet = function(){   
			console.log("Go away!");
		}

		m.rectPath = function(cx, cy, L, H){
			var rect = [[cx - L/2, cy - H/2], [cx +L/2, cy - H/2], [cx + L/2, cy + H/2], [cx - L/2, cy + H/2], [cx - L/2, cy - H/2]];
			return this.pathFromPoints(rect);
		}

		m.circlePath = function(cx, cy, r, n){
			R = [];
			theta = this.linspace(0, 2*PI, n);
			for(i=0; i<n; i++){
				R.push(r);
			}
			var x,y;
			[x,y] = this.polarToCart(R, theta);
			var xx, yy;
			[xx,yy] = this.center(x, y, cx, cy);
			return this.makeLineData(xx,yy);
		}


		m.clipImage = function(d3, svg, img, clipPath, interp, name, attr){
			if(typeof(attr) == 'undefined'){attr={};}
			if(typeof(name)=='none'){name = "temp";}
			if(typeof(interp)=='undefined'){interp="linear";}
			if(typeof(attr.x)=='undefined'){attr.x = 0;}
			if(typeof(attr.y) == 'undefined'){attr.y = 0;}
			if(typeof(attr.height) == 'undefined'){attr.heigt = svg.attr("height");}
			if(typeof(attr.width) == 'undefined'){attr.width = svg.attr("width");}


			var lineFunc = this.linefunction(d3, interp);

			svg.append("g").attr("id", name);
			pathFunction = this.linefunction(d3, interp);

			svg.append("clipPath")
				.attr("id", "clip-"+name)
				.append("path")
				.attr("d", pathFunction(clipPath));

			svg.select("#"+name).selectAll("."+name)
				.data([0])
				.enter()
				.append("image")
				.attr("xlink:href", img)
				.attr("clip-path", "url(#clip-"+name+")")
				.attr(attr);
		}

		m.clipPath = function(d3, svg, path, name, interp){
			if(typeof(name) == 'undefined'){name = "temp";}
			if(typeof(interp) == 'undefined'){interp = "linear";}
			pathFunction = this.linefunction(d3, interp);
			svg.append("clipPath")
				.attr("id", "clip-"+name)
				.append("path")
				.attr("d", pathFunction(path));
		}

		m.pathFromPoints = function(points, op){
			var xx = [];
			var yy = [];

			for(pp=0; pp<points.length; pp++){
				xx.push(points[pp][0]);
				yy.push(points[pp][1]);
			}
			return this.makeLineData(xx, yy, op)
		}

		m.polarToCart = function(r, theta){
			var x = [];
			var y = [];
			for(i=0; i<r.length; i++){
				x.push(r[i]* Math.cos(theta[i]));
				y.push(r[i]* Math.sin(theta[i]));
			}
			return [x,y]
		}

		m.makesvg = function(d3, h, w){
			if(typeof(h)=='undefined'){h=500;}
			if(typeof(w)=='undefined'){w=500;}

			return d3.select("body").append("svg")
				.attr("height", h)
				.attr("width", w);
		}

		m.makeLineData = function(xdata, ydata, op){
			if(typeof(op) == 'undefined'){
				var linedat = [];
				for(i=0; i<xdata.length; i++){
					linedat.push({x:xdata[i], y:ydata[i]});
				}
				return linedat;
			}
			else{
				var linedat = [];
				for(i=0; i<xdata.length; i++){
					linedat.push({x:xdata[i], y:ydata[i]});
				}
				var line = {p:linedat, o:op};
				return line;
			}
		}

		m.linspace = function(start, stop, size){
			var l = [];
			var step = (stop - start)/(size-1);
			for(i=0; i<(size-1); i++){
				l.push(start+i*step);
			}
			l.push(stop);
			return l;
		}

		m.area = function(points){
			if(typeof(points.p) != 'undefined'){ points = points.p;}

			var area = 0,
            aa,
            bb,
            point1,
            point2;

	        for (aa = 0, bb = points.length - 1; aa < points.length; bb=aa,aa++) {
	            point1 = points[aa];
	            point2 = points[bb];
	            area += point1.x * point2.y;
	            area -= point1.y * point2.x;
	        }
        	area /= 2;
        	return area;
		}

		m.movePath = function(path, dx, dy){
			if(typeof(path.p) != 'undefined'){
				var newpath = [];
				for(l=0; l<path.p.length; l++){
					newpath.push({x:(path.p[l].x + dx), y:(path.p[l].y + dy)});
				}

				path.p = newpath;
			}
			else{
				var newpath = [];
				for(l=0; l<path.length; l++){
					newpath.push({x:(path[l].x + dx), y:(path[l].y + dy)});
				}

				path = newpath;

			}
			return path;
		}

		m.scalePath = function(path, scalefactor){
			var cx, cy;
			[cx, cy] = this.centroidPath(path);

			if(typeof(path.p) != 'undefined'){
				var N = path.p.length;
				var ppath = path.p;
			}
			else{
				var N = path.length;
				var ppath = path;
			}
			
			var x = [], y = [];

			for(sp = 0; sp<N; sp++){
				x.push(scalefactor*(ppath[sp].x - cx) + cx);
				y.push(scalefactor*(ppath[sp].y - cy) + cy);
			}

			if(typeof(path.p) != 'undefined'){
				path.p = this.makeLineData(x,y);
				return path;
			}
			path = this.makeLineData(x, y);
			return path; 
		}

		m.rotatePath = function(path, theta, cx, cy){
			if(typeof(cx) == 'undefined' || typeof(cy) == 'undefined'){
				[cx, cy] = this.centroidPath(path);
			}
			if(typeof(path.p) != 'undefined'){
				var newpath = [];
				var xx, yy;
				for(l=0; l<path.p.length; l++){
					[xx, yy] = this.rotate(cx, cy, path.p[l].x, path.p[l].y, theta);
					newpath.push({x:xx, y:yy});
				}

				path.p = newpath;
			}			
			else{
				var newpath = [];
				var xx, yy;
				for(l=0; l<path.length; l++){
					[xx, yy] = this.rotate(cx, cy, path[l].x, path[l].y, theta);
					newpath.push({x:xx, y:yy});
				}

				path = newpath;

			}
			return path;
		}

		m.movePoints = function(points, dx, dy){
			var newpoints = [];

			for(mp=0;mp<points.length; mp++){
				newpoints.push([points[mp][0]+dx, points[mp][1]+dy]);
			}

			return newpoints;
		}

		m.plotPaths = function(d3, svg, paths, interp, name, attr){
			var ppath
			if(typeof(name) == 'undefined'){name = "temp";}
			if(typeof(interp) == 'undefined'){interp = "linear";}
			if(typeof(attr.stroke) == 'undefined'){stroke = "black";}
			if(typeof(attr["stroke-width"]) == 'undefined'){strokewidth = 1;}
			if(typeof(attr.opacity) == 'undefined'){opacity = 1;}
			if(typeof(paths[0].p) != 'undefined'){var p = true;}
			else{var p = false;}

			var lineFunc = this.linefunction(d3, interp);
			if(p){ lineFunc = function(d){return lineFunc(d.p);}}

			svg.append("g").attr("id", name);

			svg.select("#"+name).selectAll("."+name)
				.data(paths)
				.enter().append("path")
				.attr("d", lineFunc)
				.attr(attr);
		}

		m.centroidPath = function(path){
			if(typeof(path.p) != 'undefined'){ 
				var points = path.p;
			}
			else{var points = path;}
	        var cx = 0,
	            cy = 0,
	            cc,
	            dd,
	            f,
	            point1,
	            point2;

	        for (cc = 0, dd = points.length - 1; cc < points.length; dd=cc,cc++) {
	            point1 = points[cc];
	            point2 = points[dd];
	            f = point1.x * point2.y - point2.x * point1.y;
	            cx += (point1.x + point2.x) * f;
	            cy += (point1.y + point2.y) * f;
	        }

	        f = this.area(points) * 6;

	        return [cx / f, cy / f];
		}

		m.plotCircles = function(d3, svg, pathPoints, name, attr){
			if(typeof(attr) == 'undefined'){attr = {};}
			if(typeof(attr.r) == 'undefined'){attr.r = 10;}
			if(typeof(attr.fill) == 'undefined'){attr.fiill = "black";}
			if(typeof(name) == 'undefined'){name = "circ_temp";}

			svg.append("g").attr("id", name);
			svg.select("#"+name).selectAll("."+name)
				.enter().append("circle")
				.attr("cx", function(d){return d.x;})
				.attr("cy", function(d){return d.y;})
				.attr(attr);
		}

		m.roughBackground = function(w, h, dx, N){
			var lines = [];
			var n=0;
			var a = 1/(dx*dx);

			while(n < N){
				var x = [], y = [];
				var xx = Math.random()*w;
				var yy = Math.random()*h;
				var dd = Math.random()*dx;

				if(Math.random() < Math.exp(-2*dd)){
					for(i=0; i<5; i++){
						x.push(dd*Math.random()+xx);
						y.push(dd*Math.random()+yy);
					}
					var o = 0.3*Math.random();
					lines.push(this.makeLineData(x, y, o));
					n++;
				}

			}
			return lines;
		}

		m.backgroundTexture = function(svg, lineFunction, N){
			var w = svg.attr("width");
			var h = svg.attr("height");

			var lines = minimal.roughBackground(w, h, 10, N);

			svg.append("g").attr("id", "background");
			svg.append("g").attr("id", "rough");

			svg.select("#rough").selectAll(".rough")
				.data(lines)
				.enter()
				.append("path")
				.attr("d", function(d){return lineFunction(d.p);})
				.attr("line-width", 2)
				.attr("fill", "white")
				.attr("opacity", function(d){return d.o;});

			svg.select("#background").append("rect")
				.attr("x", 0)
				.attr("y", 0)
				.attr("width", "100%")
				.attr("height", "100%")
				.attr("fill", "#383838");

		}

		m.backgroundIn = function(path, size, N, w, h){
			var ppath;

			if(typeof(path.p) != 'undefined'){ppath = path.p;}
			else{ppath = path;}

			var polyX, polyY;

			[polyX, polyY] = this.pathToPoly(ppath);
			var len = polyX.length;

			var lines = [];
			var n=0;
			var a = 1/(size*size);

			while(n < N){
				var x = [], y = [];
				var xx = Math.random()*w;
				var yy = Math.random()*h;
				var dd = (1 - 2*Math.random())*size;

				if(this.inside(polyX, polyY, len, xx, yy)){
					for(i=0; i<5; i++){
						x.push(dd*Math.random()+xx);
						y.push(dd*Math.random()+yy);
					}
					var o = 0.3*Math.random();
					lines.push(this.makeLineData(x, y, o));
					n++;
				}

			}

			return lines;
		}

		m.spotFill = function(d3, svg, paths, fill, N, size){
			if(typeof(N) == 'undefined'){ N = 5000;}
			if(typeof(size) == 'undefined'){size = 3;}
			if(typeof(fill) == 'undefined'){fill = "black";}
			var w = svg.attr("width");
			var h = svg.attr("height");

			var linefunction = this.linefunction(d3);
			var specfunction = this.linefunction(d3, "basis-closed");

			var specks = this.backgroundIn(path, size, N, w, h);

			svg.append("g").attr("id", "background");
			svg.append("g").attr("id", "specks");

			svg.select("#background")
				.append("path")
				.attr("d", linefunction(path))
				.attr("fill", "black");


			svg.select("#specks").selectAll(".specks")
				.data(specks).enter()
				.append("path")
				.attr("d", function(d){return specfunction(d.p);})
				.attr("fill", "white")
				.attr("opacity", function(d){return d.o;});
			
		}

		m.spotFillAll = function(d3, svg, paths, fill, N, size){
			if(typeof(N) == 'undefined'){ N = 5000;}
			if(typeof(size) == 'undefined'){size = 3;}
			if(typeof(fill) == 'undefined'){fill = "black";}
			var w = svg.attr("width");
			var h = svg.attr("height");

			var linefunction = this.linefunction(d3);
			var specfunction = this.linefunction(d3, "basis-closed");

			var specks = this.backgroundIn(paths[0], size, N, w, h);
			var speck;

			for(pp =1; pp<paths.length; pp++){
				speck = this.backgroundIn(paths[pp], size, N, w, h);
				for(ss=0;ss<speck.length; ss++){
					specks.push(speck[ss]);
				}
			}


			svg.append("g").attr("id", "background");
			svg.append("g").attr("id", "specks");

			svg.select("#background").selectAll(".backgrouond")
				.data(paths).enter()
				.append("path")
				.attr("d", linefunction)
				.attr("fill", fill);

				console.log(speck);
			svg.select("#specks").selectAll(".specks")
				.data(specks).enter()
				.append("path")
				.attr("d", function(d){return specfunction(d.p);})
				.attr("fill", "white")
				.attr("opacity", function(d){return d.o;});
		}

		m.plotPath = function(d3, svg, path, interp, name, attr){
			var ppath
			if(typeof(name) == 'undefined'){name = "temp";}
			if(typeof(interp) == 'undefined'){interp = "linear";}
			if(typeof(attr) == 'undefined'){attr = {};}
			if(typeof(attr.stroke) == 'undefined'){attr.stroke = "black";}
			if(typeof(attr["stroke-width"]) == 'undefined'){attr["stroke-width"] = 1;}
			if(typeof(attr.fill)=='undefined'){attr.fill = "none";}
			if(typeof(attr.opacity) == 'undefined'){opacity = 1;}
			if(typeof(path.p) != 'undefined'){ppath = path.p; var p = true;}
			else{ppath = path; var p = false;}

			var lineFunc = this.linefunction(d3, interp);
			svg.append("g").attr("id", name);
			svg.select("#"+name).append("path")
				.attr("d", lineFunc(ppath))
				.attr(attr);

		}

		m.linefunction = function(d3, interpolate){
			if(typeof(interpolate) == 'undefined'){ var interpolate = 'linear';}
			return d3.svg.line().x(function(d){return d.x;}).y(function(d){return d.y;})
				.interpolate(interpolate);
		}

		m.spirals = function(cx, cy, r, p0, p1, loops){
			if(typeof(loops) == 'undefined'){ loops = 3;}
			if(typeof(p0) == 'undefined'){ p0 = 2;}
			if(typeof(p1) == 'undefined'){ p1 = 3.5;}

			var len = loops*(2*Math.PI);
			var theta = this.linspace(0, len, Math.floor(loops*1000));
			var r1 = [];
			var r2 = [];

			for(i=0; i<theta.length; i++){
				r1.push(r_1(theta[i]));
				r2.push(r_2(theta[i]));
			}

			var x1_, y1_, x2_, y2_;
			var x1, y1, x2, y2;

			[x1_, y1_] = this.polarToCart(r1, theta);
			[x2_, y2_] = this.polarToCart(r2, theta);

			[x1, y1] = this.center(x1_, y1_, cx, cy);
			[x2, y2] = this.center(x2_, y2_, cx, cy);

			var lines = this.makeLineData(x1, y1).concat(this.makeLineData(x2, y2).reverse());

			function r_1(theta){
				return r*Math.pow(theta/len, p0);
			}

			function r_2(theta){
				return r*Math.pow(theta/len, p1);
			}


			return lines;
		}

		m.randomPoints = function(w, h, N){
			var points = [];

			for(i=0; i<N; i++){
				points.push([w*Math.random(), h*Math.random()]);
			}
			return points;
		}

		m.center = function(x, y, cx, cy){
			if(typeof(cx) == 'object'){
				var h = svg.attr("height");
				var w = svg.attr("width");
			}
			else{
				var w = cx*2;
				var h = cy*2;
			}


			var newx = [],
				newy = [];
			for(i=0; i< x.length; i++){
				newx.push(x[i] + w/2);
				newy.push(y[i] + h/2);
			}
			return [newx, newy];
		}

		m.spreadFunc = function(r, R, power){
			if(typeof(power) == 'undefined'){ power = 4;}

			var rho = Math.sqrt(Math.abs(r));
			var p = Math.pow(rho/R, power);
			var rand = Math.random();
			return rand < p;
		}

		m.pointCircle = function(R, rx, ry, N, f){
			if(typeof(f) == 'undefined'){ f = this.spreadFunc;}
			if(typeof(N) == 'undefined'){ N = 1000;}

			var points = [];
			var i = 0;
			var x, y, r, o;
			while(i < N){
				x = (R+rx) - 2*R*Math.random();
				y = (R+ry) - 2*R*Math.random();
				r = (rx-x)*(rx-x)+(ry -y)*(ry - y);
				if(r < R*R && f(r, R)){
					o = Math.random()/5;
					points.push([x,y,o]);
					i++;
				}
			}
			return points;
		}

		m.plotPointCircle = function(svg, R, rx, ry, N, f){
			if(typeof(rx) == 'undefined'){ rx = svg.attr("width")/2;}
			if(typeof(ry) == 'undefined'){ ry = svg.attr("height")/2;}
			if(typeof(R) == 'undefined'){ R = svg.attr("width")/10;}

			var points = this.pointCircle(R, rx, ry, N, f);

			svg.append("g").attr("id", "dots");

			svg.select("#dots").selectAll(".dots")
			   .data(points)
			   .enter()
			   .append("circle")
			   .attr("cx", function(d){return d[0];})
			   .attr("cy", function(d){return d[1];})
			   .attr("opacity", function(d){return d[2];})
			   .attr("fill", "black")
			   .attr("r", function(d){ return Math.random();});
		}

		m.spiralInc = function(cx, cy, R){
		}

		m.rotate = function(cx, cy, x, y, angle) {
    		var radians = (Math.PI / 180) * angle,
        	cos = Math.cos(radians),
        	sin = Math.sin(radians),
        	nx = (cos * (x - cx)) + (sin * (y - cy)) + cx,
        	ny = (cos * (y - cy)) - (sin * (x - cx)) + cy;
    		return [nx, ny];
		}

		m.rotateArray = function(x, y, cx, cy, theta){
			var newx = [];
			var newy = [];
			var xx, yy;

			for(i=0; i<x.length; i++){
				[xx, yy] = this.rotate(cx, cy, x[i], y[i], theta*180/Math.PI);
				newx.push(xx);
				newy.push(yy);
			}
			return [newx, newy];
		}

		m.line = function(p0, p1, N){
			if(typeof(N)=='undefined'){ N=500;}
			var x0 = p0[0];
			var y0 = p0[1];
			var x1 = p1[0];
			var y1 = p1[1];

			function x(t){
				return x0*(1-t) + x1*t;
			}
			function y(t){
				return y0*(1-t) + y1*t;
			}
			var dt= 1/N;
			var xline = [],
				yline = [];

			for(i=0; i<N;i++){
				xline.push(x(i*dt));
				yline.push(y(i*dt));
			}
			return [xline, yline];
		}

		m.graph = function(points, connectedness){
			if(typeof(connectedness) == 'undefined'){ connectedness = 1;}

			var Nodes = [];
			var numNodes = points.length;
			var lines = [];

			for(l=0; l<numNodes; l++){
				var newNode = {};
				newNode.center = points[l];
				newNode.connections = [];
				Nodes.push(newNode);
			}

			var ck, cj;
			for(k=0; k<numNodes; k++){
				for(j = k; j<numNodes; j++){
					if(k != j){
						ck = Nodes[k].connections.length;
						cj = Nodes[j].connections.length;
						if(Math.random() < connectedness && ck < 4 && cj < 4){
							lines.push(this.line(points[k], points[j], 100));
							Nodes[k].connections.push(points[j]);
							Nodes[j].connections.push(points[k]);
						}
					}
				}
			}
			for(m=0; m<numNodes; m++){
				ck = Nodes[m].connections.length;
				while(ck<1){
					rand = Math.floor(Math.random()*numNodes);
					if(rand != m){
						lines.push(this.line(points[m], points[rand], 100));
						Nodes[m].connections.push(points[rand]);
						Nodes[rand].connections.push(points[m]);
						ck++;
					}
				}
			}

			return [points, lines];
		}

		m.inside = function(polyX, polyY, polyCorners, x, y){
			var j = polyCorners-1;
			var oddNodes = false;

			for(i=0; i<polyCorners; i++){
				if(polyY[i]<y && polyY[j]>=y
				|| polyY[j]<y && polyY[i]>=y
				&& (polyX[i]<=x || polyX[j]<=x)){
					if (polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x) {
        				oddNodes=!oddNodes; 
        			}
				}
				j=i;
			}
			return oddNodes;
		}

		m.pathToPoly = function(path){
			var ppath;
			if(typeof(path.p) != 'undefined'){ppath = path.p;}
			else{ppath = path;}

			var polyX = [];
				polyY = [];

			for(ptp = 0; ptp< ppath.length; ptp++){
				polyX.push(ppath[ptp].x);
				polyY.push(ppath[ptp].y);
			}

			return [polyX, polyY];
		}

		m.triangulate = function(p0, p1, p2, N){
			 var x0 = p0[0];
			 var y0 = p0[1];

			 var x1 = p1[0];
			 var y1 = p1[1];

			 var x2 = p2[0];
			 var y2 = p2[1];

			 var wt = Math.abs(x2 - x1);
			 var ht = Math.abs(y1 - y0);

			 var wspace = this.linspace(x1, x2, N);
			 var hspace = this.linspace(y1, y0, N);
			 var num_base = N;
			 var thislevel, x_start, x_stop;
			 var levels = [];

			 for(r=0; r<hspace.length; r++){
			 	thislevel = hspace[r];
			 	x_start = leftside(thislevel);
			 	x_stop = rightside(thislevel);
			 	levels.push(this.linspace(x_start, x_stop, num_base));
			 	num_base--;
			 }

			 var triangles = [];

			 for(k=0; k<hspace.length-1; k++){
			 	for(m=0; m<levels[k].length-1; m++){
			 		triangles.push([[levels[k][m], hspace[k]], [levels[k+1][m], hspace[k+1]], [levels[k][m+1], hspace[k]]]);
			 		if(k>0){
			 			triangles.push([[levels[k][m], hspace[k]], [levels[k-1][m+1], hspace[k-1]], [levels[k][m+1], hspace[k]]]);
			 		}
			 	}
			 }

			 var xx = [], yy=[];
			 var lines = [];

			 for(p=0; p<triangles.length; p++){
			 	xx = [];
			 	yy = [];
			 	for(r=0; r<3; r++){
			 		xx.push(triangles[p][r][0]);
			 		yy.push(triangles[p][r][1]);
			 	}
			 	lines.push([xx, yy]);
			 }

			 function leftside(y){
			 	return (x0 - x1)*y/(y0 - y1) - (x0*y1 - y0*x1)/(y0 - y1);
			 }

			 function rightside(y){
			 	return (x0 - x2)*y/(y0 - y2) - (x0*y2 - y0*x2)/(y0 - y2);
			 }

			 return lines;
		}
		
		m.triangulateToPath = function(triangles){
			var paths = [];
			var path;
			for(tt=0; tt<triangles.length; tt++){
				path = [];
				for(tp=0; tp < 3; tp++){
					path.push({x:triangles[tt][0][tp], y:triangles[tt][1][tp]});
				}
				paths.push(path);
			}
			return paths;
		}

		m.leaf = function(cx, cy, L, W, theta){
			var golden = 0.5*(1+Math.sqrt(5.0));
			theta = theta + 90;

			var x1 = cx + W/2;
			var y1 = cy + L/golden;

			var x2 = cx;
			var y2 = cy + L;

			var x3 = cx - W/2;
			var y3 = y1;

			[x1, y1] = this.rotate(cx, cy, x1, y1, theta);
			[x2, y2] = this.rotate(cx, cy, x2, y2, theta);
			[x3, y3] = this.rotate(cx, cy, x3, y3, theta);

			var points = [[cx, cy], [x1, y1], [x2, y2], [x3, y3], [cx, cy]];

			var xpoints = [cx, x1, x2, x3, cx], 
				ypoints = [cy, y1, y2, y3, cy];

			return [xpoints, ypoints];
		}
			
		return m;
	}

	if(typeof(m) == 'undefined'){
		window.m = define_minimal();
	}
	else{
		console.log("minimal already defined");
	}
})(window)