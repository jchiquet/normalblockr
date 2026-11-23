// pkgdown rewrites the links to inherited R6 methods as "<Class>.html#method-<name>",
// while roxygen2 names the anchors "method-<Class>-<name>": restore the full anchor.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll("span.pkg-link[data-topic][data-id] > a").forEach(function (a) {
    var span = a.parentElement;
    var href = a.getAttribute("href");
    var short = "#method-" + span.dataset.id;
    if (href && href.endsWith(short)) {
      a.setAttribute("href", href.slice(0, -short.length) + "#method-" + span.dataset.topic + "-" + span.dataset.id);
    }
  });
});
