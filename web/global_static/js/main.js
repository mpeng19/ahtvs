
var VOTE_CLASSES = {"undefined": "",
                    "-1": "vote-no",
                    "0": "vote-meh",
                    "1": "vote-yes",}


function draw_cards(candidate_data, vote_data) {
    var template = $('#card-template').html();
    var nav_template = $('#navcard-template').html();
    Mustache.parse(template);   // optional, speeds up future uses
    Mustache.parse(nav_template);   // optional, speeds up future uses

    var items = [];
    var context = {};
    if (candidate_data.meta.next) {
        context["next"] = candidate_data.meta.next.split("?")[1];
    }
    if (candidate_data.meta.previous) {
        context["previous"] = candidate_data.meta.previous.split("?")[1];
    }

    context["lower"] = candidate_data.meta.offset;
    context["upper"] = candidate_data.meta.offset + candidate_data.objects.length;
    context["total_count"] = candidate_data.meta.total_count;

    navcard = Mustache.to_html(nav_template, context);
    if (candidate_data.meta.next || candidate_data.meta.previous) {
        items.push(navcard);
    }

    $.each( candidate_data.objects, function(i, mol)   {
        if (typeof vote_data !== "undefined" && typeof vote_data[mol.id] !== "undefined") {
            mol["voting"] = true;
            mol["vote_pk"] = vote_data[mol.id].id;
            var rating = vote_data[mol.id].rating;
            mol["vote_rating"] = rating;
            var vote_class = VOTE_CLASSES[rating];
            var yes_class = "";
            var meh_class = "";
            var no_class = "";
            if (vote_class === "vote-yes") {yes_class = "btn-success";} else {yes_class="btn-default"}
            if (vote_class === "vote-no") {no_class = " btn-danger";} else {no_class="btn-default"}
            if (vote_class === "vote-meh") {meh_class = "btn-warning";} else {meh_class="btn-default"}

            mol["yes_class"] = yes_class;
            mol["no_class"] = no_class;
            mol["meh_class"] = meh_class;

        }
        mol["project"] = mol.project;
        mol["acceptor"] = mol.acceptor;
        mol["donor"] = mol.donor;
        mol["id"] = mol.id

        mol["weight"] = parseFloat(mol.weight).toFixed(2);
        mol["sascore"] = parseFloat(mol.sascore).toFixed(1);

        mol["splitting"] = parseFloat(mol.splitting).toFixed(3);
        mol["strength"] = parseFloat(mol.strength).toFixed(3);
        mol["absorption"] = parseFloat(mol.absorption).toFixed(2);
        mol["rate"] = parseFloat(mol.rate).toFixed(2);
        mol["homo"] = parseFloat(mol.homo).toFixed(2);
        mol["lumo"] = parseFloat(mol.lumo).toFixed(2);
        mol["inchi_key"] = mol.inchi_key;   
        mol["shortened_inchi_key"] = mol.inchi_key.split("-")[0];
        var nm = absorption_to_nm(mol.absorption);
        var rgb = Math.nmToRGB(nm);
        nm = Math.round (nm/5) * 5

        mol["color_red"] = rgb.red;
        mol["color_green"] = rgb.green;
        mol["color_blue"] = rgb.blue;
        mol["color_nm"] = nm;
        var rendered = Mustache.to_html(template, mol);

        items.push( rendered );
    });
    if (candidate_data.meta.next || candidate_data.meta.previous) {
        items.push(navcard);
    }
    $('#molecules').html(items.join( "" ));

    $(".prev-next-button").click(function(e) {
        event.preventDefault();
        var new_query = this.attributes["data-query"].value;
        load_hash(new_query);
    });

    var id_list = $.map(candidate_data.objects, function(d) {return d.id;})
    var id_string = id_list.join(",")
    var comments_by_candidate;

    $.getJSON( "/rf/comments/?format=json&candidate="+id_string, function (comments) {
        comments_by_candidate = comments.results.reduce(function(o, v, i) {
                  if (v.candidate in o) {
                    o[v.candidate].push(v)
                  } else {
                    o[v.candidate] = [v]
                  }
                  return o;
                }, {});
        $.each(comments_by_candidate, function(candidate_id, comment_group) {
            $("#"+candidate_id+" .comment-count").text(comment_group.length)
            $("#"+candidate_id+" .comment-tag").show()
        })

    });

}

function load_cards(query, ballot_id) {
//   console.log(query);
   $("#progress").show()
   $.getJSON( "/api/latest/candidate/?" + query, function( candidate_data ) {
        if (typeof ballot_id !== "undefined") {
            $.getJSON( "/api/latest/vote/?ballot=" + ballot_id, function( vote_data ) {
                votes_by_candidate = vote_data.objects.reduce(function(o, v, i) {
                      o[v.candidate.split("/")[4]] = v;
                      return o;
                    }, {});
                draw_cards(candidate_data, votes_by_candidate);
                $(".vote-button").click(function() {
                    $.post(window.location.pathname + "/vote", {
                        action: this.id
                    })
                        .done(function(data) {
                            var vote_id = data.id;
                            var vote_rating = data.rating;
                            var vote_class = VOTE_CLASSES[vote_rating];
                            var yes_class = "";
                            var meh_class = "";
                            var no_class = "";
                            if (vote_class === "vote-yes") {yes_class = "btn-success";} else {yes_class="btn-default"}
                            if (vote_class === "vote-no") {no_class = " btn-danger";} else {no_class="btn-default"}
                            if (vote_class === "vote-meh") {meh_class = "btn-warning";} else {meh_class="btn-default"}

                            $("#up-" + vote_id).removeClass().addClass("vote-button btn " + yes_class);
                            $("#meh-" + vote_id).removeClass().addClass("vote-button btn " + meh_class);
                            $("#down-" + vote_id).removeClass().addClass("vote-button btn " + no_class);

                            update_tallies(data.tallies)
                        });
                });


            });
        } else {
            draw_cards(candidate_data);
        }
        $("#progress").hide()
    });
}

function rf_load_cards(query, ballot_id) {
//   console.log(query);
   $("#progress").show()
   $.getJSON( "/rf/candidates/?" + query, function( candidate_data ) {
        if (typeof ballot_id !== "undefined") {
            $.getJSON( "/api/latest/vote/?ballot=" + ballot_id, function( vote_data ) {
                votes_by_candidate = vote_data.objects.reduce(function(o, v, i) {
                      o[v.candidate.split("/")[4]] = v;
                      return o;
                    }, {});
                rf_draw_cards(candidate_data, votes_by_candidate);
                $(".vote-button").click(function() {
                    $.post(window.location.pathname + "/vote", {
                        action: this.id
                    })
                        .done(function(data) {
                            var vote_id = data.id;
                            var vote_rating = data.rating;
                            var vote_class = VOTE_CLASSES[vote_rating];
                            var yes_class = "";
                            var meh_class = "";
                            var no_class = "";
                            if (vote_class === "vote-yes") {yes_class = "btn-success";} else {yes_class="btn-default"}
                            if (vote_class === "vote-no") {no_class = " btn-danger";} else {no_class="btn-default"}
                            if (vote_class === "vote-meh") {meh_class = "btn-warning";} else {meh_class="btn-default"}

                            $("#up-" + vote_id).removeClass().addClass("vote-button btn " + yes_class);
                            $("#meh-" + vote_id).removeClass().addClass("vote-button btn " + meh_class);
                            $("#down-" + vote_id).removeClass().addClass("vote-button btn " + no_class);

                            update_tallies(data.tallies)
                        });
                });


            });
        } else {
            rf_draw_cards(candidate_data);
        }
        $("#progress").hide()
    });
}

function get_url_vars_after_hash()
{
    var vars = {}, hash;
    var hashes = window.location.hash.substring(1).split('&');
    for(var i = 0; i < hashes.length; i++)
    {
        hash = hashes[i].split('=');
        vars[hash[0]] = hash[1];
    }
    return vars;
}

function rf_draw_cards(candidate_data, vote_data) {
    var template = $('#card-template').html();
    var nav_template = $('#navcard-template').html();
    Mustache.parse(template);   // optional, speeds up future uses
    Mustache.parse(nav_template);   // optional, speeds up future uses

    var items = [];
    var context = {};
    if (candidate_data.next) {
        context["next"] = candidate_data.next.split("?")[1];
    }
    if (candidate_data.previous) {
        context["previous"] = candidate_data.previous.split("?")[1];
    }

    url_vars = get_url_vars_after_hash()
    if (url_vars.offset) {
        context["lower"] = parseInt(url_vars.offset, 10);
    }
    else {
        context["lower"] = 0;
    }


    method = url_vars.method || "b3lyp_tddft_631gs_rpa_s0_geom";

    context["upper"] = Math.min(context["lower"] + candidate_data.results.length, candidate_data.count);

    context["total_count"] = candidate_data.count;

    navcard = Mustache.to_html(nav_template, context);
    if (candidate_data.next || candidate_data.previous) {
        items.push(navcard);
    }


    $.each( candidate_data.results, function(i, mol)   {
        if (typeof vote_data !== "undefined" && typeof vote_data[mol.id] !== "undefined") {
            mol["voting"] = true;
            mol["vote_pk"] = vote_data[mol.id].id;
            var rating = vote_data[mol.id].rating;
            mol["vote_rating"] = rating;
            var vote_class = VOTE_CLASSES[rating];
            var yes_class = "";
            var meh_class = "";
            var no_class = "";
            if (vote_class === "vote-yes") {yes_class = "btn-success";} else {yes_class="btn-default"}
            if (vote_class === "vote-no") {no_class = " btn-danger";} else {no_class="btn-default"}
            if (vote_class === "vote-meh") {meh_class = "btn-warning";} else {meh_class="btn-default"}

            mol["yes_class"] = yes_class;
            mol["no_class"] = no_class;
            mol["meh_class"] = meh_class;

        }
        mol["project"] = mol.project;
        mol["acceptor"] = mol.acceptor;
        mol["donor"] = mol.donor;
        mol["id"] = mol.id
        mol["inchi_key"] = mol.inchi_key;   
        mol["shortened_inchi_key"] = mol.inchi_key.split("-")[0];
        mol["weight"] = parseFloat(mol.weight).toFixed(2);
        mol["sascore"] = parseFloat(mol.sascore).toFixed(1);

        props = {}
        mol.properties.filter(function(p) {return p.method == "b3lyp_tddft_631gs_rpa_s0_geom"})
            .forEach(function (p) {props[p.name] = p.value;})

        mol["splitting"] = parseFloat(props.splitting).toFixed(3);
        mol["strength"] = parseFloat(props.strength).toFixed(3);
        mol["absorption"] = parseFloat(props.absorption).toFixed(2);
        mol["homo"] = parseFloat(props.homo).toFixed(2);
        mol["lumo"] = parseFloat(props.lumo).toFixed(2);
        mol["rate"] = parseFloat(props.rate).toFixed(2);


        var nm = absorption_to_nm(mol.absorption);
        var rgb = Math.nmToRGB(nm);
        nm = Math.round (nm/5) * 5

        mol["color_red"] = rgb.red;
        mol["color_green"] = rgb.green;
        mol["color_blue"] = rgb.blue;
        mol["color_nm"] = nm;
        var rendered = Mustache.to_html(template, mol);

        items.push( rendered );
    });
    if (candidate_data.next || candidate_data.previous) {
        items.push(navcard);
    }
    $('#molecules').html(items.join( "" ));

    $(".prev-next-button").click(function(e) {
        event.preventDefault();
        var new_query = this.attributes["data-query"].value;
        var hashes = new_query.split('&');
        add_query("offset=", false)
        add_query("limit=", false)
        for(var i = 0; i < hashes.length; i++)
        {
            add_query(hashes[i], false);
        }
        var hash = window.location.hash.substring(1);
        load_hash(hash);
    });

    var id_list = $.map(candidate_data.results, function(d) {return d.id;})
    var id_string = id_list.join(",")
    var comments_by_candidate;

    $.getJSON( "/rf/comments/?format=json&candidate="+id_string, function (comments) {
        comments_by_candidate = comments.results.reduce(function(o, v, i) {
                  if (v.candidate in o) {
                    o[v.candidate].push(v)
                  } else {
                    o[v.candidate] = [v]
                  }
                  return o;
                }, {});
        $.each(comments_by_candidate, function(candidate_id, comment_group) {
            $("#"+candidate_id+" .comment-count").text(comment_group.length)
            $("#"+candidate_id+" .comment-tag").show()
        })

    });

}


 $.fn.buildHarvey = function () {
        
        var size    = parseInt(this.data('size')),
            radius  = size / 2.0,
            value   = parseInt(this.data('value'));
        
        // pie all the things!
        if (this.length > 1){
            this.each( function() {
                $(this).buildHarvey();
            });
            return this;
        }
        
        // is you trying to break things?
        if (isNaN(value)) {
            return this;
        }
        
        // set the size of this
        this.css({
            height: size,
            width: size
        }).addClass('pie-sliced');
        
        if (value >= 95) {
            this.addClass('highlight');
        }
        // make value behave
        value = Math.min(Math.max(value, 0), 100);

        // make me some svg
        this.pie = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        
        // if value is 100 or higher, just use a circle
        if (value >= 100) {
            this.pie.slice = document.createElementNS("http://www.w3.org/2000/svg", "circle");
            this.pie.slice.setAttribute('r', radius);
            this.pie.slice.setAttribute('cx', radius);
            this.pie.slice.setAttribute('cy', radius);
            
        } else {
            this.pie.slice = document.createElementNS("http://www.w3.org/2000/svg", "path");
            
            //calculate x,y coordinates of the point on the circle to draw the arc to. 
            var x = Math.cos((2 * Math.PI)/(100/value));
            var y = Math.sin((2 * Math.PI)/(100/value));
            
            //should the arc go the long way round?
            var longArc = (value <= 50) ? 0 : 1;
            
            //d is a string that describes the path of the slice.
            var d = "M" + radius + "," + radius + " L" + radius + "," + 0 + ", A" + radius + "," + radius + " 0 " + longArc + ",1 " + (radius + y*radius) + "," + (radius - x*radius) + " z";       
            this.pie.slice.setAttribute('d', d);
        }
        
        //add the slice to the pie.
        $(this.pie.slice).appendTo(this.pie);
        
        // add the pie to this
        $(this.pie).appendTo(this);
        
        return this;
    };

function load_redox_cards(query) {
    $("#progress").show()
    var template = $('#card-template').html();
    var nav_template = $('#navcard-template').html();
    Mustache.parse(template);   // optional, speeds up future uses
    Mustache.parse(nav_template);   // optional, speeds up future uses

    $.getJSON( "/api/latest/redoxpair/?" + query, function( data ) {
        var items = [];
        var context = {};
        if (data.meta.next) {
            context["next"] = data.meta.next.split("?")[1];
        }
        if (data.meta.previous) {
            context["previous"] = data.meta.previous.split("?")[1];
        }

        context["lower"] = data.meta.offset;
        context["upper"] = data.meta.offset + data.objects.length;
        context["total_count"] = data.meta.total_count;

        navcard = Mustache.to_html(nav_template, context);
        if (data.meta.next || data.meta.previous) {
            items.push(navcard);
        }
        $.each( data.objects, function(i, redox)   {
            redox["id"] = redox.id;
            redox["weight"] = parseFloat(redox.reduced.weight).toFixed(2);
            redox["sascore"] = parseFloat(redox.reduced.sascore).toFixed(1);
            var rp = parseFloat(redox.redox_potential);
            redox["redox_potential"] = rp.toFixed(2);
            redox['high_redox_potential_percent'] = Math.round(Math.max(0, Math.min(Math.round((1 - Math.abs(1 - rp) * 4) * 100), 100)));
            redox['low_redox_potential_percent'] = Math.round(Math.max(0, Math.min(Math.round((1 - Math.abs(0 - rp) * 4) * 100), 100)));
            redox["solvation"] = parseFloat(redox.reduced.water_solvation_energy).toFixed(2);
            redox['solvation_percent'] = Math.round((1 - Math.abs(1 - rp) * 4) * 100);

            redox["homo"] = parseFloat(redox.reduced.homo).toFixed(2);
            redox["lumo"] = parseFloat(redox.reduced.lumo).toFixed(2);
            redox["red_shortened_inchi_key"] = redox.reduced.inchi_key.split("-")[0];
            redox["red_inchi_key"] = redox.reduced.inchi_key;
            redox["ox_shortened_inchi_key"] = redox.oxidized.inchi_key.split("-")[0];
            redox["ox_inchi_key"] = redox.oxidized.inchi_key;
            redox["red_cas"] = redox.reduced.cas;
            redox["ox_cas"] = redox.oxidized.cas;
            redox["red_sa"] = parseFloat(redox.reduced.sascore).toFixed(1);

            var logK = parseFloat(redox.log_hyd_constant)
            if (logK > 0.00001 || logK < -0.00001) {
                redox["log_hyd_constant"] = logK.toFixed(2);
                redox['log_hyd_constant_percent'] = Math.round(Math.max(0, Math.min((3 - logK) * 33.33, 100)))
            } else {
                redox["log_hyd_constant"] = 'NA'
            }

            var mike = parseFloat(redox.michael_hydration_energy);
            if (mike > 0.00001 || mike < -0.00001) {
                redox["michael_hydration_energy"] = mike.toFixed(2);
                redox['michael_hydration_percent'] = Math.round(Math.max(0, Math.min(mike * 200 + 100, 100)))
            } else {
                redox["michael_hydration_energy"] = 'NA'
            }

            if (redox["michael_hydration_energy"] == "NaN") redox["michael_hydration_energy"] = "";
            if (redox["red_sa"] == "NaN") redox["red_sa"] = "";

            var rendered = Mustache.to_html(template, redox);

            items.push( rendered );
        });
        if (data.meta.next || data.meta.previous) {
            items.push(navcard);
        }
        $('#molecules').html(items.join( "" ));

        $(".prev-next-button").click(function(e) {
            event.preventDefault();
            var new_query = this.attributes["data-query"].value;
            load_hash(new_query);
        });
        $('.harvey').buildHarvey();
        $("#progress").hide()
    });
}



function update_tallies(tallies) {
    percent_up = 100.0 * tallies.up / tallies.total
    percent_down = 100.0 * tallies.down / tallies.total
    percent_meh = 100.0 * tallies.meh / tallies.total
    percent_remaining = 100.0 * tallies.remaining / tallies.total

    $("#progress-remaining").width(percent_remaining + "%")
    $("#progress-up").width(percent_up + "%")
    $("#progress-down").width(percent_down + "%")
    $("#progress-meh").width(percent_meh + "%")

    $("#total-up").text(tallies["up"])
    $("#total-down").text(tallies["down"])
    $("#total-meh").text(tallies["meh"])
    $("#total-remaining").text(tallies["remaining"])
    //$("#total-done").text(tallies["total"] - tallies["remaining"])
    $("#total").text(tallies["total"])
}
