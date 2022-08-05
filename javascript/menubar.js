/** Show or hide title bar menu */
function showMenu(){
    var menu = document.getElementById("menu");
    var btn  = document.getElementById("menu-btn");
    if(menu.className === "my-menu"){
        menu.className += " responsive";
        btn.innerHTML   = "&times;"
    }
    else {
        menu.className = "my-menu";
        btn.innerHTML  = "☰";
    }
}