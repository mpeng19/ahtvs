var all_tags = {};

var margin = {
    top: 50,
    right: 20,
    bottom: 10,
    left: 80
}

var twoPi = 2 * Math.PI,
    progress = 0,
    total = 1308573, // must be hard-coded if server doesn't report Content-Length
    formatPercent = d3.format(".0%");

var window_height = 800;

var _global_layout;
var xhr;

var colorscale = d3.scale.quantize()
                    .domain([4,-4]).range(colorbrewer.RdBu[4])

var oled_layout = {"x_label": "Oscillator Strength",
                   "x": function(d) {return d.strength;},
                   "y_label": "Splitting (eV)",
                   "y": function(d) {return d.splitting;},
                   "r": function(d) {return d.rate;},
                   "color":  function(d) {rgb = Math.nmToRGB(absorption_to_nm(d.absorption));
                                                return d3.rgb(rgb.red, rgb.green, rgb.blue);
                                                },
                   "card": get_oled_card,

                   "minx": 0,
                   "maxx": 2,
                   "miny": 0,
                   "maxy": 2,
                   "minr": 0, 
                   "maxr": 0.5
                  }


var redox_layout = {"x_label": "Redox Potential (eV)",
                    "x": function(d) {return d.redox_potential;},
                    "y_label": "Solvation Energy (eV)",
                     "y": function (d) {return d.water_solvation_energy},
                    "r": function(d) { if (d.log_hyd_constant === 0) return 0.2; else return 0.3;},
                    "color":  function(d) {val = d.log_hyd_constant;
                                              if (val === 0) return "rgba(100,100,100,0.4)";
                                              return colorscale(val);
                                            },
                     "card": get_redox_card,
                     "minx": -1,
                     "maxx": 2,
                     "miny": -3.5,
                     "maxy": 0,
                     "minr": 0, 
                     "maxr": 0.5
                     }

function load_url(url, layout) {
    $("#progress").show()  // show the spinner
    // if there is already an open request running, abort and start a new one.
    // this is for when we update filters in the interactive plots, and want to stop the
    // current request
    if (typeof xhr !== "undefined") {
        xhr.abort();
        console.log("abort")
    }
    xhr = d3.json(url)
        // success handler
        .on("load", function(input_data) { 
            data = input_data.objects;
            update_bubbles(data, layout);
            $("#progress").hide()  // hide the spinner
        })
        // error handler
        .on("error", function(error) {  
            console.log("failure!", error); 
        })
        .get();

}

function load_rf_url(url, layout) {
    $("#progress").show()  // show the spinner
    // if there is already an open request running, abort and start a new one.
    // this is for when we update filters in the interactive plots, and want to stop the
    // current request
    if (typeof xhr !== "undefined") {
        xhr.abort();
        console.log("abort")
    }
    xhr = d3.json(url)
        // success handler
        .on("load", function(input_data) { 
            data = input_data.results;
            data.forEach(function(mol) {
                mol.properties.filter(function(p) {return p.method == "b3lyp_tddft_631gs_rpa_s0_geom"})
                    .forEach(function (p) {mol[p.name] = parseFloat(p.value).toFixed(3);})
            });

            update_bubbles(data, layout);
            $("#progress").hide()  // hide the spinner
        })
        // error handler
        .on("error", function(error) {  
            console.log("failure!", error); 
        })
        .get();

}


function load_bubbles(query, layout) {
    load_rf_url("/rf/candidates/?format=json&limit=5000&"+query, layout);
}

function load_redox(query, layout) {
    load_url("/api/latest/miniredoxpair/?format=json&"+query, layout);
}

function get_oled_card(div, d) {
    div.append("a")
        .attr("href", "/detail/"+d.project+"/"+d.inchi_key)
        .append("img")
         .attr("src", image_url(d.inchi_key))
         .attr("class", "card-main-img");
    var table = div.append("table")
                .attr("class", "table table-striped")
                .append("tbody");
    d3.json("/rf/candidates/"+d.id, function(d) {
        d.properties.filter(function(p) {return p.method == "b3lyp_tddft_631gs_rpa_s0_geom"})
                    .forEach(function (p) {d[p.name] = parseFloat(p.value).toFixed(3);})
        name_data_list = [
            ["Nicknames", d.nicknames],
            ["Splitting", d.splitting + " (eV)"],
            ["Strength", d.strength],
            ["Absorption", d.absorption + " (eV)"],
            ["Rate", d.rate + " (1/us)"],
            ["HOMO", d.homo + " (eV)"],
            ["LUMO", d.lumo + " (eV)"],
            ["SA Score", d.sascore],
            ["Key", d.inchi_key.split("-")[0]]
        ]

        table.selectAll("tr")
            .data(name_data_list)
            .enter()
            .append("tr")
            .each(function(r) {
                var tr = d3.select(this)
                tr.append("td").text(r[0])
                tr.append("td").text(r[1])
            });
        var lastrow = table.append("tr")
        lastrow.append("td").text("Color Est.: ")
        var colortd = lastrow.append("td")
        var nm = absorption_to_nm(d.absorption);
        nm = Math.round (nm/5) * 5
        colortd.append("span").text("  ~" + nm.toString() + "nm ")
        colortd.append("svg")
            .attr("width", "16")
            .attr("height", "14")
            .append("circle")
            .attr("r", "7")
            .attr("cx", "7")
            .attr("cy", "7")
            .style("fill", function() {
                rgb = Math.nmToRGB(nm);
                return d3.rgb(rgb.red, rgb.green, rgb.blue);
            })
    })
}

function get_redox_card(div, pair) {
    d3.json(pair.resource_uri, function(d) {
        name_data_list = [
            ["Nicknames", d.reduced.nicknames],
            ["Redox Potential", parseFloat(d.redox_potential).toFixed(2)],
            ["Solvation E", parseFloat(d.water_solvation_energy).toFixed(2)],
            ["Log K Hyd", parseFloat(d.log_hyd_constant).toFixed(2)],
            ["Key", d.reduced.inchi_key.split("-")[0]]
        ]

        div.append("img")
            .attr("src", '/api/latest/candidate/'+d.reduced.id+'?format=svg')
            .attr("class", "card-main-img");
        div.append("img")
            .attr("src", '/api/latest/candidate/'+d.oxidized.id+'?format=svg')
            .attr("class", "card-main-img");
        var table = div.append("table")
                    .attr("class", "table table-striped")
                    .append("tbody");
        table.selectAll("tr")
            .data(name_data_list)
            .enter()
            .append("tr")
            .each(function(r) {
                var tr = d3.select(this)
                tr.append("td").text(r[0])
                tr.append("td").text(r[1])
            });
    })
}


function init_oled_bubbles(query) {
    init_bubbles(oled_layout);
    if (typeof query !== "undefined") {
        load_bubbles(query, oled_layout);
    }
}

function init_redox_bubbles(query) {
    init_bubbles(redox_layout);
    if (typeof query !== "undefined") {
        load_redox(query, redox_layout);
    }
}


function init_bubbles(layout) {
    if (typeof(svg) === "undefined") {

        window_width = parseInt(d3.select('#chart').style('width'), 10) - margin.left - margin.right;

        window_xlen = Math.min(layout.maxx, (layout.maxy / window_height) * window_width);

        d3.select(window).on('resize', get_resize_func_for_layout(layout));


        xscale = d3.scale.linear()
            .domain([layout.minx, layout.maxx])
            .range([0, window_width]);

        yscale = d3.scale.linear()
            .domain([layout.miny, layout.maxy])
            .range([0, window_height]);

        rscale = d3.scale.linear()
            .domain([layout.minr, layout.maxr])
            .range([1, 4])

        xAxis = d3.svg.axis()
            .scale(xscale)
            .orient("top")
            .tickFormat(d3.format(".2f"))
            .ticks(5)

        yAxis = d3.svg.axis()
            .scale(yscale)
            .orient("left")
            .tickFormat(d3.format(".2f"))
            .ticks(5)


        zoom = d3.behavior.zoom()
            .scaleExtent([1, 400])
            .x(xscale)
            .y(yscale)
            .on("zoom", get_zoom_func_for_layout(layout));

        svg = d3.select("#chart").append("svg")
            .attr("width", window_width + margin.left + margin.right + 20)
            .attr("height", window_height + margin.top + margin.bottom)

        //.style("margin-left", margin.left + "px")

        svg.append("text")
            .attr("class", "x axlabel")
            .attr("text-anchor", "middle")
            .attr("x", window_width / 2)
            .attr("y", "1em")
            .text(layout["x_label"]);

        svg.append("text")
            .attr("class", "y axlabel")
            .attr("text-anchor", "middle")
            .attr("x", -window_height / 2)
            .attr("y", "1em")
            .attr("transform", "rotate(-90)")
            .text(layout["y_label"]);

        chart = svg
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
            .attr("class", "chart")

        imgarea = chart.append("g")
            .attr("class", "background")
            .attr("width", window_width)
            .attr("height", window_height)
            .attr('clip-path', 'url(#clip-grid)')
            .call(zoom);

        // add background to catch mouse events for zooming
        bg = imgarea.append("rect")
            .attr("class", "chart-background")
            .attr("height", window_height)
            .attr("width", window_width)

        clipGrid = chart.append("clipPath")
            .attr('id', 'clip-grid')
            .append('rect')
            .attr('x', 0)
            .attr('y', 0)
            .attr('width', window_width)
            .attr('height', window_height);

        chart.append("g")
            .attr("class", "x axis")
            .call(xAxis);

        chart.append("g")
            .attr("class", "y axis")
            .call(yAxis);

        border = svg.append("rect")
            .attr("class", "chart-border")
            .attr("x", margin.left)
            .attr("y", margin.top)
            .attr("height", window_height)
            .attr("width", window_width - 1)


    }

}



function get_resize_func_for_layout(layout) {
    var my_zoom = get_zoom_func_for_layout(layout);
    function resize() {
        // update width
        var window_width = parseInt(d3.select('#chart').style('width'), 10) - margin.left - margin.right;
        var yd = yscale.domain();
        var window_zoomed_domain_width = Math.min(layout.maxx - layout.minx, ((yd[1] - yd[0]) / window_height) * window_width);
        var window_full_domain_width = Math.min(layout.maxx - layout.minx, ((layout.maxy - layout.miny)/ window_height) * window_width);

        svg.attr("width", window_width + margin.left + margin.right);
        // reset x range

        clipGrid.attr("width", window_width); // update clipping

        var xd0 = xscale.domain()[0]

        xscale.domain([layout.minx, layout.minx + window_full_domain_width])
            .range([0, window_width])

        var old_scale = zoom.scale()
        var old_trans = zoom.translate()
        zoom.x(xscale);

        zoom.scale(old_scale)
        zoom.translate(old_trans)

        xscale.domain([xd0, xd0 + window_zoomed_domain_width])

        limit_domains(layout);

        chart.select(".x.axis").call(xAxis);

        svg.selectAll(".row_line line").attr("x2", window_width)
        svg.select(".x.axlabel").attr("x", window_width / 2)

        bg.attr("width", window_width);
        border.attr("width", window_width);
        my_zoom();

    } // end resize
    return resize
}

function limit_domains(layout) {
    var xd = xscale.domain();
    var yd = yscale.domain();
    var t = zoom.translate()
    var tx = t[0];
    var ty = t[1];
    var hit_limit = false;

    if (xd[0] < layout.minx) {
        xscale.domain([layout.minx, xd[1] - xd[0]]);
        tx = xscale(layout.minx);
        hit_limit = true;
    }

    if (xd[1] > layout.maxx) {
        xscale.domain([xd[0] - (xd[1] - layout.maxx), layout.maxx]);
        tx += xscale(xd[1]) - xscale(layout.maxx);
        hit_limit = true;
    }

    if (yd[0] < layout.miny) {
        yscale.domain([layout.miny, yd[1] - yd[0]]);
        ty = yscale(layout.miny);
        hit_limit = true;
    }

    if (yd[1] > layout.maxy) {
        yscale.domain([yd[0] - (yd[1] - layout.maxy), layout.maxy]);
        ty += yscale(yd[1]) - yscale(layout.maxy);
        hit_limit = true;
    }


    zoom.translate([tx, ty]);

    return hit_limit;
} // end limit_domains


function get_zoom_func_for_layout(layout) {
    function zoomed() {
        limit_domains(layout);

        chart.select(".x.axis").call(xAxis);
        chart.select(".y.axis").call(yAxis);

        chart.selectAll(".dot")
            .attr("cx", function(d) {
                return xscale(layout["x"](d));
            })
            .attr("cy", function(d) {
                return yscale(layout["y"](d));
            })
            .attr("r", function(d) {
                return (((zoom.scale() - 1) * 0.05) + 1) * rscale(layout["r"](d));
            })
    }
    return zoomed
}


function get_parameter_by_name(name) {
    name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
    var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
        results = regex.exec(location.search);
    return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

function image_url(inchi_key) {
    return '/inchi_img/' + inchi_key + ".svg"
    //return '/api/latest/molecule/?inchi_key=' + inchi_key + '&format=svg'
}


function update_bubbles(data, layout) {
    // sort by strength to put smallest circles on top
    data.sort(function(a, b) {
        return (b.rate > a.rate) ? 1 : -1;

    })

    imgarea.selectAll(".dot")
        .remove()

    imgarea.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("r", function(d) {
            if (d.splitting < 0) {
                return 0;
            }
            if (d.rate > 0) {
                return rscale(d.rate);
            } else {
                return 0;
            }
        })
        .attr("cx", function(d) {
            return xscale(layout["x"](d));
        })
        .attr("cy", function(d) {
            return yscale(layout["y"](d));
        })
        .style("fill", function(d) {
            if (d.gray) {
                return "#ccc";
            } else {
                return layout["color"](d)
            }
        })
        .on("click", function(d) {
            d3.select(".selected").classed("selected", false)
            d3.select(this).classed("selected", true);

            $("#selected-mol-card").empty()
            var details = d3.select("#selected-mol-card")

            layout["card"](details, d);


        })


    if (data.length > 0) {

        max_x = d3.max(data, function(d) {
            return parseFloat(layout["x"](d), 10);
        });

        max_y= d3.max(data, function(d) {
            return parseFloat(layout["y"](d), 10);
        });

        zoom.scale(d3.max([1, d3.min([layout.maxx / max_y, layout.maxy / max_x])* 0.95]));
        zoom.event(imgarea.transition().duration(1000));
    }





} // end filter function
