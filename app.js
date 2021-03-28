import { Main } from './src/Main.js'

function main() {
 
const container = document.querySelector('#container');    
const instance = new Main(container);
instance.render();

}

main();
