
function add_node_tooltips() {
    tippy('[data-tippy-content]', {
        allowHTML: true,
        placement: 'bottom',
        theme: 'translucent',
        interactive: true,
    });
}

document.addEventListener("DOMContentLoaded", function (event) {
    add_node_tooltips();
});